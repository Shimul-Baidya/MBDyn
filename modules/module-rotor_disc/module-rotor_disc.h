// Author: Matteo Daniele <matteo.daniele@polimi.it>
// tail rotor module for mbdyn
#ifndef MODULE_ROTOR_DISC_H
#define MODULE_ROTOR_DISC_H

extern bool RotorDiscSet(void);

#include "elem.h"
#include "dataman.h"
#include "userelem.h"
#include "strforce.h"
#include "tpldrive_impl.h"
#include "strforce_impl.h"

/*
FROM PROUTY "HELICOPTER PERFORMANCE, STABILITY AND OCONTROL, CHAPTER 3"
V1 IS CALLED CONSTANT MOMENTUM INDUCED VELOCITY: 
GOOD APPROXIMATION OF INDUCED VELOCITY AT SPEED > 30KT
*/

class RotorDisc : virtual public Elem, public UserDefinedElem
{
private:

    const StructNode* pHubNode;
    // position of the force wrt hub reference node
    Vec3 HubNodeArm;
    // the disc rotor force vector (Thrust oriented along the z of the local rf)
    Vec3 OutputThrust;

    // orientation matrix of thrust vector
    Mat3x3 RThrustOrientation;

    // control input (rad)
    const DriveCaller* pXColl;
    // air density driver (kg/m3)
    const DriveCaller* pRho;
    //rotor angular speed driver (rad/s)
    const DriveCaller* pOmega;

    // partial derivatives initialization
    doublereal dMu[5];
    doublereal dLambda[5];

    // partial derivatives of thrust components for Jacobian assembly
    doublereal dT0[5];
    doublereal dTTheta[5];
    doublereal dTLambda[5];

    // elements of thrust that will be updated at each time step
    doublereal T0;
    doublereal TTheta;
    doublereal TLambda;

    // thrust value at each time step
    doublereal Thrust;
    // induced drag value at each time step
    doublereal DragInduced;
    // induced power at each time step
    doublereal PowerInduced;
    // elements of thrust for jacobian assembly
    doublereal dThrust[5];
    
    // partial derivatives of advance ratio
    void dMuCalc(void);
    // partial derivatives of inflow ratio
    void dLambdaCalc(void);

    // partial derivatives of thrust components for jacobian assembly
    void dT0Calc(void);
    void dTThetaCalc(void);
    void dTLambdaCalc(void);

    // thrust components 
    void T0Calc(void);
    void TThetaCalc(void);
    void TLambdaCalc(void);
    // thrust and induced power (output as private data?)
    void ThrustCalc(void);

    // jacobians of thrust
    void dTCalc(void);

    // state-dependent variable calculation
    void updateStatesDeps(void);

    
    
    // dimensional data of the tail rotor
    doublereal RotorRadius;           // tail rotor radius [m]
    doublereal DiscArea;           // tail rotor disc area [m^2]
    doublereal RotorSolidity;       // tail rotor solidity [-]
    doublereal ClAlpha ;     // lift slope curve [1/rad]
    doublereal BladeTwist;    // twist angle at tip [rad]

    // Cl stall angle [rad] 
    doublereal AOAStallMin;        
    doublereal AOAStallMax;        
    
    // main rotor data for v1hover
    doublereal hubs_distance;           // distance between mr and tr hub [m]
    doublereal mr_nominal_power_shp;    // main rotor nominal power, OGE [hp]
    doublereal mr_nominal_omega;        // main rotor nominal angular speed [rad/s]

    // shp2W
    const doublereal sHP2W = 745.7; // conversion factor, horsepower to watt

    // part of v1h depending only on the constant factors:
    // v1hPart = sqrt(Th/2A), where Th is the required rotor 
    // thrust in hover to compensate for nominal rotor torque
    doublereal Th;
    doublereal v1hPart;

    // airspeed in body frame acting on tail rotor center node (x pointing towards nose, z pointing outwards)
    Vec3 VTrHub;
    // tail rotor angular speed
    doublereal RotorOmega;
    // air density
    doublereal rho;
    // STATE DEPENDENT VARIABLES
    doublereal u;
    doublereal v;
    doublereal w;
    doublereal Vtot;
    doublereal Vtot2;
    // tip speed
    doublereal Vtip;
    doublereal Vtip2;
    // advance ratio
    doublereal mu;
    doublereal mu2;
    doublereal mu4;
    // v1h: induced velocity in hover
    doublereal v1h;
    // CONSTANT MOMENTUM INDUCED VELOCITY
    doublereal a_v1; // sqrt((vtot/2)^2+v1h^4)
    doublereal V1;
    // inflow ratio
    doublereal lambda;

    // pedal input [rad]
    doublereal thetaColl;
    // pedal angle limits [rad]
    doublereal thetaCollMin; 
    doublereal thetaCollMax; 

    Vec3 F;
    Vec3 M;

public:

    // constructor inherits FollowerForceClass constructor
    RotorDisc( unsigned int uL, const DofOwner *pDO,
                    DataManager* pDM, MBDynParser& HP);
    virtual ~RotorDisc();
    
    // update inputs - check saturation limits
    void inputSaturation();
    // compute thrust
    void computeRotorThrust();
    // method to build partial derivatives and jacobians
    void assemblyJacobian();



    ////////////////////////////////////////////////////////////////////////
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
    virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
    //contributo al file di restart
    std::ostream& Restart(std::ostream& out) const;
    // assemblaggio yacobiano
    VariableSubMatrixHandler& AssJac(   VariableSubMatrixHandler& WorkMat, 
                                        doublereal dCoef,
	                                    const VectorHandler& /* XCurr */ ,
	                                    const VectorHandler& /* XPrimeCurr */ );
    // assemblaggio residuo
    SubVectorHandler& AssRes(   SubVectorHandler& WorkVec,
	                            doublereal /* dCoef */ ,
	                            const VectorHandler& /* XCurr */ ,
	                            const VectorHandler& /* XPrimeCurr */ );
    unsigned int iGetNumPrivData(void) const;
    unsigned int iGetPrivDataIdx(const char* s) const;
    doublereal dGetPrivData(unsigned int i) const;
    int iGetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);
	virtual unsigned int iGetInitialNumDof(void) const;
    // output
    virtual void Output(OutputHandler& OH) const;
    /* Contributo allo jacobiano durante l'assemblaggio iniziale */
    VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,
                                            const VectorHandler& /* XCurr */ );
    /* Contributo al residuo durante l'assemblaggio iniziale */
    SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
                                    const VectorHandler& /* XCurr */ );

};

#endif // ! MODULE_ROTOR_DISC_H