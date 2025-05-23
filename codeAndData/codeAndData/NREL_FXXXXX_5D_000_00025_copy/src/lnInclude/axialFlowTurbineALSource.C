/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author(s)
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of turbinesFoam, which is based on OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "axialFlowTurbineALSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "unitConversion.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(axialFlowTurbineALSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        axialFlowTurbineALSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::axialFlowTurbineALSource::createCoordinateSystem()
{
    // Make sure axis is a unit vector
    axis_ /= mag(axis_);

    // Free stream direction is a unit vector
    freeStreamDirection_ = freeStreamVelocity_ / mag(freeStreamVelocity_);

    // Radial direction is vertical direction
    verticalDirection_ /= mag(verticalDirection_);
    radialDirection_ = verticalDirection_;

    pitchingAxis_ = - freeStreamDirection_ ^ verticalDirection_; // Added by Clem 0202

    // Calculate initial azimuthal or tangential direction
    azimuthalDirection_ = axis_ ^ radialDirection_;             // Modified by Clem 0202
    azimuthalDirection_ /= mag(azimuthalDirection_);
}


void Foam::fv::axialFlowTurbineALSource::createBlades()
{
    int nBlades = nBlades_;
    blades_.setSize(nBlades);
    int nElements;
    List<List<scalar> > elementData;
    word modelType = "actuatorLineSource";
    List<scalar> frontalAreas(nBlades); // frontal area from each blade

    for (int i = 0; i < nBlades_; i++)
    {
        word bladeName = bladeNames_[i];
        // Create dictionary items for this blade
        dictionary bladeSubDict = bladesDict_.subDict(bladeName);
        bladeSubDict.lookup("nElements") >> nElements;
        bladeSubDict.lookup("elementData") >> elementData;
        scalar azimuthalOffset = bladeSubDict.lookupOrDefault
        (
            "azimuthalOffset",
            0.0
        );

        bladeSubDict.add("freeStreamVelocity", freeStreamVelocity_);
        bladeSubDict.add("fieldNames", coeffs_.lookup("fieldNames"));
        bladeSubDict.add("profileData", profileData_);

        // Disable individual lifting line end effects model if rotor-level
        // end effects model is active
        if
        (
            bladeSubDict.found("endEffects")
            and endEffectsActive_
            and endEffectsModel_ != "liftingLine"
        )
        {
            bladeSubDict.set("endEffects", false);
        }
        else if (endEffectsModel_ == "liftingLine" and endEffectsActive_)
        {
            bladeSubDict.add("endEffects", true);
        }

        if (debug)
        {
            Info<< "Creating actuator line blade " << bladeName << endl;
            Info<< "Blade has " << nElements << " elements" << endl;
            Info<< "Element data:" << endl;
            Info<< elementData << endl << endl;
        }

        // Convert element data into actuator line element geometry
        label nGeomPoints = elementData.size();
        List<List<List<scalar> > > elementGeometry(nGeomPoints);
        List<vector> initialVelocities(nGeomPoints, vector::zero);
        // Frontal area for this blade
        scalar frontalArea = 0.0;
        scalar maxRadius = 0.0;
        forAll(elementData, j)
        {
            // Read AFTAL dict element data
            scalar axialDistance = elementData[j][0];
            scalar radius = elementData[j][1];
            scalar azimuthDegrees = elementData[j][2] + azimuthalOffset;
            scalar azimuthRadians = degToRad(azimuthDegrees);
            scalar chordLength = elementData[j][3];
            scalar chordMount = elementData[j][4];
            scalar pitch = elementData[j][5];

            // Find max radius for calculating frontal area
            if (radius > maxRadius)
            {
                maxRadius = radius;
            }

            // Set sizes for actuatorLineSource elementGeometry lists
            elementGeometry[j].setSize(6);
            elementGeometry[j][0].setSize(3);
            elementGeometry[j][1].setSize(3);
            elementGeometry[j][2].setSize(1);
            elementGeometry[j][3].setSize(3);
            elementGeometry[j][4].setSize(1);
            elementGeometry[j][5].setSize(1);

            // Create geometry point for AL source at origin
            vector point = origin_;
            // Move point along axial direction
            point += axialDistance*axis_;
            // Move along radial direction
            point += radius*radialDirection_;
            // Move along chord according to chordMount
            scalar chordDisplacement = (chordMount - 0.25)*chordLength;
            point -= chordDisplacement*azimuthalDirection_;
            // Set initial velocity of quarter chord
            scalar radiusCorr = sqrt(magSqr((chordMount - 0.25)*chordLength)
                                     + magSqr(radius));
            vector initialVelocity = azimuthalDirection_*omega_*radiusCorr;
            scalar velAngle = atan2(((chordMount - 0.25)*chordLength), radius);
            rotateVector(initialVelocity, vector::zero, axis_, velAngle);
            initialVelocities[j] = initialVelocity;
            // Rotate point and initial velocity according to azimuth value
            rotateVector(point, origin_, axis_, azimuthRadians);
            rotateVector
            (
                initialVelocities[j],
                vector::zero,
                axis_,
                azimuthRadians
            );

            // Set point coordinates for AL source
            elementGeometry[j][0][0] = point.x(); // x location of geom point
            elementGeometry[j][0][1] = point.y(); // y location of geom point
            elementGeometry[j][0][2] = point.z(); // z location of geom point

            // Set chord reference direction
            vector chordDirection = azimuthalDirection_;

            // Set span directions for AL source
            // Blades start oriented vertically
            // Use planform normal to figure out span direction
            vector planformNormal = freeStreamDirection_;
            vector spanDirection = chordDirection ^ planformNormal;
            spanDirection /= mag(spanDirection);

            // Rotate span and chord directions according to azimuth
            rotateVector(spanDirection, vector::zero, axis_, azimuthRadians);
            elementGeometry[j][1][0] = spanDirection.x();
            elementGeometry[j][1][1] = spanDirection.y();
            elementGeometry[j][1][2] = spanDirection.z();
            rotateVector(chordDirection, vector::zero, axis_, azimuthRadians);
            elementGeometry[j][3][0] = chordDirection.x();
            elementGeometry[j][3][1] = chordDirection.y();
            elementGeometry[j][3][2] = chordDirection.z();

            // Set chord length
            elementGeometry[j][2][0] = chordLength;

            // Set chord mount
            elementGeometry[j][4][0] = chordMount;

            // Set element pitch or twist
            elementGeometry[j][5][0] = pitch;
        }

        // Add frontal area to list
        frontalArea = mathematical::pi*magSqr(maxRadius);
        frontalAreas[i] = frontalArea;

        if (debug)
        {
            Info<< "Converted element geometry:" << endl << elementGeometry
                << endl;
            Info<< "Frontal area from " << bladeName << ": " << frontalArea
                << endl;
        }

        bladeSubDict.add("elementGeometry", elementGeometry);
        bladeSubDict.add("initialVelocities", initialVelocities);
        bladeSubDict.add("dynamicStall", dynamicStallDict_);
        bladeSubDict.add
        (
            "velocitySampleRadius",
            coeffs_.lookupOrDefault("velocitySampleRadius", 0.0)
        );
        bladeSubDict.add
        (
            "nVelocitySamples",
            coeffs_.lookupOrDefault("nVelocitySamples", 20)
        );
        bladeSubDict.add("selectionMode", coeffs_.lookup("selectionMode"));
        bladeSubDict.add("cellSet", coeffs_.lookup("cellSet"));

        // Do not write force from individual actuator line unless specified
        bladeSubDict.lookupOrAddDefault("writeForceField", false);


        bladeSubDict.add("sixDofDataAL", sixDofData_);  // Added by Clem 0123
        bladeSubDict.add("epsilonFactor", epsilonFactor_);  // Added by Clem 0131

        dictionary dict;
        dict.add("actuatorLineSourceCoeffs", bladeSubDict);
        dict.add("type", "actuatorLineSource");
        dict.add("active", dict_.lookup("active"));

        actuatorLineSource* blade = new actuatorLineSource
        (
            name_ + "." + bladeName,
            modelType,
            dict,
            mesh_
        );

        blades_.set(i, blade);
    }

    // Frontal area is calculated using defined rotorRadius rather than
    // detected from elementData
    frontalArea_ = mathematical::pi*magSqr(rotorRadius_);
    Info<< "Frontal area of " << name_ << ": " << frontalArea_ << endl;
}


void Foam::fv::axialFlowTurbineALSource::createHub()
{
    int nElements;
    List<List<scalar> > elementData;
    dictionary hubSubDict = hubDict_;

    hubDict_.lookup("nElements") >> nElements;
    hubDict_.lookup("elementData") >> elementData;

    // Convert element data into actuator line element geometry
    label nGeomPoints = elementData.size();
    List<List<List<scalar> > > elementGeometry(nGeomPoints);
    List<vector> initialVelocities(nGeomPoints, vector::zero);

    forAll(elementData, j)
    {
        // Read hub element data
        scalar axialDistance = elementData[j][0];
        scalar height = elementData[j][1];
        scalar diameter = elementData[j][2];

        // Set sizes for actuatorLineSource elementGeometry lists
        elementGeometry[j].setSize(6);
        elementGeometry[j][0].setSize(3);
        elementGeometry[j][1].setSize(3);
        elementGeometry[j][2].setSize(1);
        elementGeometry[j][3].setSize(3);
        elementGeometry[j][4].setSize(1);
        elementGeometry[j][5].setSize(1);

        // Create geometry point for AL source at origin
        vector point = origin_;
        // Move along axis
        point += axialDistance*axis_;
        // Move along vertical direction
        point += height*verticalDirection_;

        elementGeometry[j][0][0] = point.x(); // x location of geom point
        elementGeometry[j][0][1] = point.y(); // y location of geom point
        elementGeometry[j][0][2] = point.z(); // z location of geom point

        // Set span directions
        elementGeometry[j][1][0] = verticalDirection_.x();
        elementGeometry[j][1][1] = verticalDirection_.y();
        elementGeometry[j][1][2] = verticalDirection_.z();

        // Set chord length
        elementGeometry[j][2][0] = diameter;

        // Set chord reference direction
        elementGeometry[j][3][0] = freeStreamDirection_.x();
        elementGeometry[j][3][1] = freeStreamDirection_.y();
        elementGeometry[j][3][2] = freeStreamDirection_.z();

        // Set chord mount
        elementGeometry[j][4][0] = 0.25;

        // Set pitch
        elementGeometry[j][5][0] = 0.0;
    }

    hubSubDict.add("elementGeometry", elementGeometry);
    hubSubDict.add("initialVelocities", initialVelocities);
    hubSubDict.add("fieldNames", coeffs_.lookup("fieldNames"));
    hubSubDict.add("profileData", profileData_);
    hubSubDict.add("freeStreamVelocity", freeStreamVelocity_);
    hubSubDict.add("selectionMode", coeffs_.lookup("selectionMode"));
    hubSubDict.add("cellSet", coeffs_.lookup("cellSet"));

    // Do not write force from individual actuator line unless specified
    hubSubDict.lookupOrAddDefault("writeForceField", false);

    hubSubDict.add("sixDofDataAL", sixDofData_);  // Added by Clem 0125
    hubSubDict.add("epsilonFactor", epsilonFactor_);  // Added by Clem 0131

    dictionary dict;
    dict.add("actuatorLineSourceCoeffs", hubSubDict);
    dict.add("type", "actuatorLineSource");
    dict.add("active", dict_.lookup("active"));

    actuatorLineSource* hub = new actuatorLineSource
    (
        name_ + ".hub",
        "actuatorLineSource",
        dict,
        mesh_
    );

    hub_.set(hub);
}


void Foam::fv::axialFlowTurbineALSource::createTower()
{
    vector towerAxis = verticalDirection_;
    List<List<scalar> > elementData;
    dictionary towerSubDict = towerDict_;

    towerDict_.lookup("elementData") >> elementData;

    // Convert element data into actuator line element geometry
    label nGeomPoints = elementData.size();
    List<List<List<scalar> > > elementGeometry(nGeomPoints);
    List<vector> initialVelocities(nGeomPoints, vector::zero);

    forAll(elementData, j)
    {
        // Read tower element data
        scalar axialDistance = elementData[j][0];
        scalar height = elementData[j][1];
        scalar diameter = elementData[j][2];

        // Set sizes for actuatorLineSource elementGeometry lists
        elementGeometry[j].setSize(6);
        elementGeometry[j][0].setSize(3);
        elementGeometry[j][1].setSize(3);
        elementGeometry[j][2].setSize(1);
        elementGeometry[j][3].setSize(3);
        elementGeometry[j][4].setSize(1);
        elementGeometry[j][5].setSize(1);

        // Create geometry point for AL source at origin
        vector point = origin_;
        // Move along turbine axis
        point += axialDistance*axis_;
        // Move along tower axis according to height
        point += height*towerAxis;

        elementGeometry[j][0][0] = point.x(); // x location of geom point
        elementGeometry[j][0][1] = point.y(); // y location of geom point
        elementGeometry[j][0][2] = point.z(); // z location of geom point

        // Set span directions
        elementGeometry[j][1][0] = towerAxis.x(); // x component of span
        elementGeometry[j][1][1] = towerAxis.y(); // y component of span
        elementGeometry[j][1][2] = towerAxis.z(); // z component of span

        // Set chord length
        elementGeometry[j][2][0] = diameter;

        // Set chord reference direction
        elementGeometry[j][3][0] = freeStreamDirection_.x();
        elementGeometry[j][3][1] = freeStreamDirection_.y();
        elementGeometry[j][3][2] = freeStreamDirection_.z();

        // Set chord mount
        elementGeometry[j][4][0] = 0.25;

        // Set pitch
        elementGeometry[j][5][0] = 0.0;
    }

    towerSubDict.add("elementGeometry", elementGeometry);
    towerSubDict.add("initialVelocities", initialVelocities);
    towerSubDict.add("fieldNames", coeffs_.lookup("fieldNames"));
    towerSubDict.add("profileData", profileData_);
    towerSubDict.add("freeStreamVelocity", freeStreamVelocity_);
    towerSubDict.add("selectionMode", coeffs_.lookup("selectionMode"));
    towerSubDict.add("cellSet", coeffs_.lookup("cellSet"));

    // Do not write force from individual actuator line unless specified
    towerSubDict.lookupOrAddDefault("writeForceField", false);

    towerSubDict.add("sixDofDataAL", sixDofData_);  // Added by Clem 0125
    towerSubDict.add("epsilonFactor", epsilonFactor_);  // Added by Clem 0131

    dictionary dict;
    dict.add("actuatorLineSourceCoeffs", towerSubDict);
    dict.add("type", "actuatorLineSource");
    dict.add("active", dict_.lookup("active"));

    actuatorLineSource* tower = new actuatorLineSource
    (
        name_ + ".tower",
        "actuatorLineSource",
        dict,
        mesh_
    );

    tower_.set(tower);
}


void Foam::fv::axialFlowTurbineALSource::createNacelle()
{
    // Do nothing
}


void Foam::fv::axialFlowTurbineALSource::calcEndEffects()
{
    if (debug)
    {
        Info<< "Calculating end effects for " << name_ << endl;
    }

    // Info<< "inside end effect 1" << name_ << endl;  // Added By Clem 0202

    // Calculate rotor-level end effects correction
    scalar pi = Foam::constant::mathematical::pi;
    // Info<< "inside end effect 2" << name_ << endl;  // Added By Clem 0202
    forAll(blades_, i)
    {
        forAll(blades_[i].elements(), j)
        {
            // Info<< "inside end effect 3" << blades_[i].elements()[j].name() << endl;  // Added By Clem 0202
            scalar rootDist = blades_[i].elements()[j].rootDistance();
            vector relVel = blades_[i].elements()[j].relativeVelocity();
            vector elementVel = blades_[i].elements()[j].velocity();
            if (debug)
            {
                Info<< "    rootDist: " << rootDist << endl;
                Info<< "    relVel: " << relVel << endl;
            }
            // Calculate angle between rotor plane and relative velocity
            scalar phi = pi/2.0;
            if (mag(relVel) > VSMALL)
            {
                vector elementVelDir = elementVel / mag(elementVel);
                scalar relVelOpElementVel = -elementVelDir & relVel;
                vector rotorPlaneDir = axis_;                     // Modified by Clem 0202 (was freeStreamDirection_)
                scalar relVelRotorPlane = rotorPlaneDir & relVel;
                // Note: Does not take yaw into account
                phi = atan2(relVelRotorPlane, relVelOpElementVel);
                phi = max(0.01, phi);                             // Added by Clem 0203 to avoid numerical error (the float cannot handle super small number).
            }
            // Info<< "phi of " << blades_[i].elements()[j].name() << " is: " << phi << endl;  // Added By Clem 0202
            if (debug)
            {
                scalar phiDeg = Foam::radToDeg(phi);
                Info<< "    phi (degrees): " << phiDeg << endl;
            }
            // Calculate end effect factor for this element
            scalar f = 1.0;
            dictionary endEffectsCoeffs = endEffectsDict_.subOrEmptyDict
            (
                endEffectsModel_ + "Coeffs"
            );
            if (endEffectsModel_ == "Glauert")
            {
                if (endEffectsCoeffs.lookupOrDefault("tipEffects", true))
                {
                    scalar acosArg = Foam::exp
                    (
                        -nBlades_/2.0*(1.0/rootDist - 1)/sin(phi)
                    );
                    // Info<< "acosArg of " << blades_[i].elements()[j].name() << " is: " << acosArg << endl << endl;  // Added By Clem 0202
                    f = 2.0/pi*acos(min(1.0, acosArg));
                }
                if (endEffectsCoeffs.lookupOrDefault("rootEffects", false))
                {
                    scalar tipDist = 1.0 - rootDist;
                    scalar acosArg = Foam::exp
                    (
                        -nBlades_/2.0*(1.0/tipDist - 1)/sin(phi)
                    );
                    f *= 2.0/pi*acos(min(1.0, acosArg));
                }
            }
            else if (endEffectsModel_ == "Shen")
            {
                scalar c1;
                endEffectsCoeffs.lookup("c1") >> c1;
                scalar c2;
                endEffectsCoeffs.lookup("c2") >> c2;
                scalar g = Foam::exp(-c1*(nBlades_*tipSpeedRatio_ - c2)) + 0.1;
                if (endEffectsCoeffs.lookupOrDefault("tipEffects", true))
                {
                    scalar acosArg = Foam::exp
                    (
                        -g*nBlades_/2.0*(1.0/rootDist - 1)/sin(phi)
                    );
                    f = 2.0/pi*acos(min(1.0, acosArg));
                }
                if (endEffectsCoeffs.lookupOrDefault("rootEffects", false))
                {
                    scalar tipDist = 1.0 - rootDist;
                    scalar acosArg = Foam::exp
                    (
                        -g*nBlades_/2.0*(1.0/tipDist - 1)/sin(phi)
                    );
                    f *= 2.0/pi*acos(min(1.0, acosArg));
                }
            }
            // Info<< "inside end effect 5" << blades_[i].elements()[j].name() << endl;  // Added By Clem 0202
            if (debug)
            {
                Info<< "    f: " << f << endl;
            }
            blades_[i].elements()[j].setEndEffectFactor(f);
        }
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::axialFlowTurbineALSource::axialFlowTurbineALSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    turbineALSource(name, modelType, dict, mesh),
    hasHub_(false),
    hasTower_(false),
    hasNacelle_(false),
    verticalDirection_
    (
        coeffs_.lookupOrDefault("verticalDirection", vector(0, 0, 1))
    )
{
    read(dict);   //Also sets some parameters, Added by Clem 0123
    createCoordinateSystem();
    createBlades();
    Info<< endl << "hasTower_: " << hasTower_  << endl << endl; // Added by Clem 0123 (checking setting of hasTower_)
    if (hasHub_)
    {
        createHub();
    }
    if (hasTower_)
    {
        createTower();
    }
    if (hasNacelle_)
    {
        createNacelle();
    }
    createOutputFile();

    // Rotate turbine to azimuthalOffset if necessary
    scalar azimuthalOffset = coeffs_.lookupOrDefault("azimuthalOffset", 0.0);
    rotate(degToRad(azimuthalOffset));

    // Yaw turbine to a static value if specified
    scalar yawAngle = coeffs_.lookupOrDefault("yawAngle", 0.0);
    yaw(degToRad(yawAngle));


    //Start added by Clem 0202
    // Pitch turbine to a static value if specified, 
    scalar pitchAngleRotor = coeffs_.lookupOrDefault("tiltAngle", 0.0);
    if (pitchAngleRotor == 0.0)
    {
        Info << "No pitch" << endl << endl;
    }
    else
    {
        pitchRotor(degToRad(pitchAngleRotor));
    }
    //End added by Clem 0202


    if (true)
    {
        Info<< "axialFlowTurbineALSource created at time = " << time_.value()
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::axialFlowTurbineALSource::~axialFlowTurbineALSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::axialFlowTurbineALSource::rotate(scalar radians)
{
    if (debug)
    {
        Info<< "Rotating " << name_ << " " << radians << " radians"
            << endl << endl;
    }

    forAll(blades_, i)
    {
        blades_[i].rotate(origin_, axis_, radians);
        blades_[i].setSpeed(origin_, axis_, omega_);
    }

    if (hasHub_)
    {
        hub_->rotate(origin_, axis_, radians);
        hub_->setSpeed(origin_, axis_, omega_);
    }
}


// Added by Clem 0202 NOTE: should be adjusted before implementing together with pitchRotor()
void Foam::fv::axialFlowTurbineALSource::yaw(scalar radians)
{
    if (debug)
    {
        Info<< "Yawing " << name_ << " " << radians << " radians"
            << endl << endl;
    }

    vector referencePoint = origin_ - origin_;

    // First, rotate axis
    rotateVector(axis_, referencePoint, verticalDirection_, radians);

    forAll(blades_, i)
    {
        blades_[i].rotate(origin_, verticalDirection_, radians);
    }

    if (hasHub_)
    {
        hub_->rotate(origin_, verticalDirection_, radians);
    }
}



// Start added by Clem 0202 NOTE: should be adjusted before implementing together with yaw()
void Foam::fv::axialFlowTurbineALSource::pitchRotor(scalar radians)
{

    bool debugClemPitch = true;

    if (debugClemPitch)
    {
        Info<< "Pitching " << name_ << " " << radians << " radians"
            << endl << endl;
        Info<< "axis_ before pitching: " << axis_ << endl << endl;
        Info<< "pitchingAxis_ before pitching: " << pitchingAxis_ << endl << endl;
    }

    vector referencePoint = origin_ - origin_;

    // First, rotate axis
    rotateVector(axis_, referencePoint, pitchingAxis_, radians);

    forAll(blades_, i)
    {
        blades_[i].rotate(origin_, pitchingAxis_, radians);
    }

    if (hasHub_)
    {
        hub_->rotate(origin_, pitchingAxis_, radians);
    }

    if (debugClemPitch)
    {
        Info<< "axis_ after pitching: " << axis_ << endl << endl;
    }

}
// End added by Clem 0202


//- Calculate the area-averaged axial induction factor Added by Clem 0214 0218
List<scalar> Foam::fv::axialFlowTurbineALSource::areaAvgInductionBEMCTCP()
{
	
	// area-averaged axial induction factor, BEM CT, BEM CP
	List<scalar> inductionBEMCTCP(nBlades_ * 3);
	
	scalar pi = Foam::constant::mathematical::pi;
	
	forAll(blades_, i)
    {
		
		scalar areaAvgInductionFac = 0;
		scalar areaAvgBEMCT = 0;
		scalar areaAvgBEMCP = 0;
		
	    scalar numericalArea = 0;
	    
		forAll(blades_[i].elements(), j)
		{
			vector bladeEleToOrig = blades_[i].elements()[j].position() - origin_;    // Get the postion vector between sctuato point and rotor center, Clem 0311
			scalar radialSpanPosition = bladeEleToOrig & blades_[i].elements()[j].spanDirection();
			scalar annulusArea = 2 * pi * radialSpanPosition * blades_[i].elements()[j].spanLength();
			numericalArea += annulusArea;
			
			//scalar axialVelocityEle = blades_[i].elements()[j].inflowVelocityOrig() & freeStreamDirection_;
			// Clem thinks that axis_ is more resonable since the assumption of BEM
			scalar axialVelocityEle = blades_[i].elements()[j].inflowVelocityOrig() & axis_;
			scalar aEle = 1 - axialVelocityEle / mag(freeStreamVelocity_);
			// axialIndcutionFactorEle
			
			scalar endEffectFactorEle = blades_[i].elements()[j].endEffectFactor();
			
			scalar BEMCTEle = 0;
			scalar BEMCPEle = 0;

			if (aEle > 1/3)
			{
				BEMCTEle = 4 * aEle * ( 1 - aEle ) * endEffectFactorEle;
				BEMCPEle = 4 * aEle * ( 1 - aEle ) * endEffectFactorEle  * ( 1 - aEle );
			}
			else
			{
				BEMCTEle = 4 * aEle * ( 1 - 0.25 * (5 - 3 * aEle) * aEle ) * endEffectFactorEle;
				BEMCPEle = 4 * aEle * ( 1 - 0.25 * (5 - 3 * aEle) * aEle ) * endEffectFactorEle * ( 1 - aEle );
			}
			
			areaAvgInductionFac += annulusArea * aEle;
			areaAvgBEMCT += annulusArea * BEMCTEle;
			areaAvgBEMCP += annulusArea * BEMCPEle;
		}
		
		inductionBEMCTCP[3*i] = areaAvgInductionFac / numericalArea;
		inductionBEMCTCP[3*i + 1] = areaAvgBEMCT / numericalArea;
		inductionBEMCTCP[3*i + 2] = areaAvgBEMCP / numericalArea;

    }
      
    return inductionBEMCTCP;
}
// End added by Clem 0214 0218




void Foam::fv::axialFlowTurbineALSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    Info<< endl; // Added by Clem 0128
    Info<< "Show something when this is called 1" << endl << endl; // Added by Clem 0127

    // Rotate the turbine if time value has changed
    if (time_.value() != lastRotationTime_)
    {
        rotate();
    }
    Info<< "Show something when finishing turb rotate" << endl << endl; // added by Clem 0202

    // Zero out force vector and field, note that each level has its forcefield variables Clem 0311
    forceField_ *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);
    force_ *= 0;

    // Check dimensions of force field and correct if necessary
    if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
    {
        forceField_.dimensions().reset(eqn.dimensions()/dimVolume);
    }

    // Create local moment vector
    vector moment(vector::zero);

    // Info<< "Before end effect" << endl; // Added by Clem 0127 // added by Clem 0202

    if (endEffectsActive_ and endEffectsModel_ != "liftingLine")
    {
        // Calculate end effects based on current velocity field
        calcEndEffects();
    }

    // Info<< "After end effect" << endl << endl; // Added by Clem 0127 // added by Clem 0202

    //Daeling with surge
    //Start Added by Clem 0127
    scalar currentTime = mesh_.time().value();
    scalar deltaT = mesh_.time().deltaT().value();
    
    // Commented out by Clem 0127
    //if (currentTime == deltaT)
    //{
    //    originOrig_ = origin_;
    //}
    
    //End Added by Clem 0127

    originPrevious_ = origin_; // Added by Clem 0127

    //Start Added by Clem 0127
    translationRotationVectors motionsixDofComplex = interpolateSplineXY//missios [C] : if timestep is within the given times in
    (
        currentTime, //the current time 
        sixDofTimeArray_, //given times in the file
        sixDofDisComplexArray_ //All the six dof displacements
    );

    sixDofDisWT_ = motionsixDofComplex[0];

    Info<< endl;
    Info<< "sixDofDisWT_ of " << name_ << " :" << sixDofDisWT_ << endl;

    sixDofVelWT_ = vector((originOrig_ + sixDofDisWT_ - originPrevious_))/deltaT; //0127
    Info<< "sixDofVelWT_ of " << name_ << " :" << sixDofVelWT_ << endl << endl; //0127

    Info<< "origin_ not changed for " << name_ << " : " << origin_ << endl;
    origin_ = originPrevious_ + sixDofVelWT_ * deltaT;
    Info<< "origin_ changed for " << name_ << " : " << origin_ << endl << endl << endl;
    //End Added by Clem 0127





    // Add source for blade actuator lines
    forAll(blades_, i)
    {
        blades_[i].addSup(eqn, fieldI, sixDofDisWT_, sixDofVelWT_);  //Modified by Clem 0127
        forceField_ += blades_[i].forceField();
        force_ += blades_[i].force();
        bladeMoments_[i] = blades_[i].moment(origin_);
        moment += bladeMoments_[i];
    }

    if (hasHub_)
    {
        // Add source for hub actuator line
        hub_->addSup(eqn, fieldI, sixDofDisWT_, sixDofVelWT_);  //Modified by Clem 0127
        forceField_ += hub_->forceField();
        force_ += hub_->force();
        moment += hub_->moment(origin_);
    }

    if (hasTower_)
    {
        // Add source for tower actuator line
        tower_->addSup(eqn, fieldI, sixDofDisWT_, sixDofVelWT_);   //Modified by Clem 0127
        forceField_ += tower_->forceField();
        if (includeTowerDrag_)
        {
            force_ += tower_->force();
        }
    }

    if (hasNacelle_)
    {
        // Add source for tower actuator line
        nacelle_->addSup(eqn, fieldI, sixDofDisWT_, sixDofVelWT_);   //Modified by Clem 0127
        forceField_ += nacelle_->forceField();
        if (includeNacelleDrag_)
        {
            force_ += nacelle_->force();
        }
    }

    // Torque is the projection of the moment from all blades on the axis
    torque_ = moment & axis_;

    torqueCoefficient_ = torque_/(0.5*frontalArea_*rotorRadius_
                       * magSqr(freeStreamVelocity_));
    powerCoefficient_ = torqueCoefficient_*tipSpeedRatio_;
    dragCoefficient_ = force_ & freeStreamDirection_
                     / (0.5*frontalArea_*magSqr(freeStreamVelocity_));

    // Print performance to terminal
    printPerf();

    // Write performance data -- note this will write multiples if there are
    // multiple PIMPLE loops
    if (Pstream::master())
    {
        writePerf();
        
        // Start Added by Clem 0214 0218
        List<scalar> inductionBEMCTCP = areaAvgInductionBEMCTCP();
        
        forAll(blades_, i)
        {
				*outputFile_<< "," << inductionBEMCTCP[3*i];
				*outputFile_<< "," << inductionBEMCTCP[3*i+1];
				*outputFile_<< "," << inductionBEMCTCP[3*i+2];
		}
		*outputFile_<< endl;
		// End Added by Clem 0214 0218
    }
}

// compressible not used Clem 0311
void Foam::fv::axialFlowTurbineALSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
     // Added by Clem 0127
    Info<< "Show something when this is called, this version do not support compressible flow" << endl << endl;
    
    // Rotate the turbine if time value has changed
    if (time_.value() != lastRotationTime_)
    {
        rotate();
    }

    // Zero out force vector and field
    forceField_ *= dimensionedScalar("zero", forceField_.dimensions(), 0.0);
    force_ *= 0;

    // Check dimensions of force field and correct if necessary
    if (forceField_.dimensions() != eqn.dimensions()/dimVolume)
    {
        forceField_.dimensions().reset(eqn.dimensions()/dimVolume);
    }

    // Create local moment vector
    vector moment(vector::zero);

    if (endEffectsActive_ and endEffectsModel_ != "liftingLine")
    {
        // Calculate end effects based on current velocity field
        calcEndEffects();
    }


    //Daeling with surge
    //Start Added by Clem 0127
    scalar currentTime = mesh_.time().value();
    scalar deltaT = mesh_.time().deltaT().value();
    if (currentTime == deltaT)
    {
        originPrevious_ = origin_;
        originOrig_ = origin_;
    }
    //End Added by Clem 0127


    //Start Added by Clem 0127
    List<Tuple2<scalar, translationRotationVectors>> data_sixDof;
    sixDofData_.readEntry("data", data_sixDof);
    translationRotationVectors motionsixDofComplex = data_sixDof[0].second();
    sixDofDisWT_ = motionsixDofComplex[0];
    Info<< "sixDofDisWT_ turbine: " << sixDofDisWT_ << endl << endl;

    Info<< "origin_ not changed turbine: " << origin_ << endl << endl;
    origin_ = originOrig_ + sixDofDisWT_;
    Info<< "origin_ changed turbine: " << origin_ << endl << endl;

    sixDofVelWT_ = vector((origin_ - originPrevious_))/deltaT; //0127
    Info<< "sixDofVelWT_ turbine: " << sixDofVelWT_ << endl << endl; //0127
    //End Added by Clem 0127

    // Add source for blade actuator lines
    forAll(blades_, i)
    {
        blades_[i].addSup(rho, eqn, fieldI);
        forceField_ += blades_[i].forceField();
        force_ += blades_[i].force();
        bladeMoments_[i] = blades_[i].moment(origin_);
        moment += bladeMoments_[i];
    }

    if (hasHub_)
    {
        // Add source for hub actuator line
        hub_->addSup(rho, eqn, fieldI);
        forceField_ += hub_->forceField();
        force_ += hub_->force();
        moment += hub_->moment(origin_);
    }

    if (hasTower_)
    {
        // Add source for tower actuator line
        tower_->addSup(rho, eqn, fieldI);
        forceField_ += tower_->forceField();
        if (includeTowerDrag_)
        {
            force_ += tower_->force();
        }
    }

    if (hasNacelle_)
    {
        // Add source for tower actuator line
        nacelle_->addSup(rho, eqn, fieldI);
        forceField_ += nacelle_->forceField();
        if (includeNacelleDrag_)
        {
            force_ += nacelle_->force();
        }
    }

    // Torque is the projection of the moment from all blades on the axis
    torque_ = moment & axis_;

    scalar rhoRef;
    coeffs_.lookup("rhoRef") >> rhoRef;
    torqueCoefficient_ = torque_/(0.5*rhoRef*frontalArea_*rotorRadius_
                       * magSqr(freeStreamVelocity_));
    powerCoefficient_ = torqueCoefficient_*tipSpeedRatio_;
    dragCoefficient_ = force_ & freeStreamDirection_
                     / (0.5*rhoRef*frontalArea_*magSqr(freeStreamVelocity_));

    // Print performance to terminal
    printPerf();

    // Write performance data -- note this will write multiples if there are
    // multiple PIMPLE loops
    if (Pstream::master())
    {
        writePerf();
    }

    originPrevious_ = origin_; // Added by Clem 0127
}


void Foam::fv::axialFlowTurbineALSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    // Rotate the turbine if time value has changed
    if (time_.value() != lastRotationTime_)
    {
        rotate();
    }

    Info<< "Show something when this is called" << endl << endl; // Added by Clem 0127

    if (endEffectsActive_ and endEffectsModel_ != "liftingLine")
    {
        // Calculate end effects based on current velocity field
        calcEndEffects();
    }

    // Add scalar source term from blades
    forAll(blades_, i)
    {
        blades_[i].addSup(eqn, fieldI);
    }

    if (hasHub_)
    {
        // Add source for hub actuator line
        hub_->addSup(eqn, fieldI);
    }

    if (hasTower_)
    {
        // Add source for tower actuator line
        tower_->addSup(eqn, fieldI);
    }

    if (hasNacelle_)
    {
        // Add source for nacelle actuator line
        nacelle_->addSup(eqn, fieldI);
    }
}


void Foam::fv::axialFlowTurbineALSource::printCoeffs() const
{
    Info<< "Number of blades: " << nBlades_ << endl;
}


bool Foam::fv::axialFlowTurbineALSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        turbineALSource::read(dict);

        // Get hub information
        hubDict_ = coeffs_.subOrEmptyDict("hub");
        if (hubDict_.keys().size() > 0)
        {
            hasHub_ = true;
        }

        // Get tower information
        towerDict_ = coeffs_.subOrEmptyDict("tower");
        if (towerDict_.keys().size() > 0)
        {
            hasTower_ = true;
        }
        includeTowerDrag_ = towerDict_.lookupOrDefault
        (
            "includeInTotalDrag",
            false
        );

        // Get nacelle information
        nacelleDict_ = coeffs_.subOrEmptyDict("nacelle");
        if (nacelleDict_.keys().size() > 0)
        {
            hasNacelle_ = true;
        }
        includeNacelleDrag_ = nacelleDict_.lookupOrDefault
        (
            "includeInTotalDrag",
            false
        );

        // Read end effects subdictionary
        endEffectsDict_ = coeffs_.subOrEmptyDict("endEffects");
        endEffectsDict_.lookup("active") >> endEffectsActive_;
        endEffectsDict_.lookup("endEffectsModel") >> endEffectsModel_;

        if (debug)
        {
            Info<< "Debugging on" << endl;
            Info<< "Axial-flow turbine properties:" << endl;
            printCoeffs();
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
