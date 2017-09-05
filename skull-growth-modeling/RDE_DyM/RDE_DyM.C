/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    RDE_DyM

Description
    Reaction-Diffusion-Strain solver with dynamic mesh technique for cranial vault growth.

Author
    Chanyoung Lee

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "constitutiveModel.H"
#include "solidInterface.H"
#include "volPointInterpolation.H"
#include "pointPatchInterpolation.H"
#include "primitivePatchInterpolation.H"
#include "pointFields.H"
#include "twoDPointCorrector.H"
#include "leastSquaresVolPointInterpolation.H"
#include "processorFvPatchFields.H"
#include "transformGeometricField.H"
#include "symmetryPolyPatch.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
# include "setRootCase.H"
# include "createTime.H"
//# include "createMesh.H"
#   include "createDynamicFvMesh.H"
# include "createFields.H"
# include "readDivDSigmaExpMethod.H"
# include "readDivDSigmaLargeStrainExpMethod.H"
# include "readMoveMeshMethod.H"
# include "createSolidInterfaceNonLin.H"
# include "findGlobalFaceZones.H"

// Read check frequency
    label checkFrequency = 1;
    args.optionReadIfPresent("checkFrequency", checkFrequency);
    
//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info << "\nStarting time loop\n" << endl;
  
  cellVolume.internalField() = mesh.V();

  for (runTime++; !runTime.end(); runTime++)
    {
      Info<< "Time = " << runTime.timeName() << nl << endl;
      


#     include "readSolidMechanicsControls.H"

cellVolume.internalField() = mesh.V();
dV = fvc::ddt(cellVolume);
rho = (rho_b-rho_m)*(pow(o,m)/(pow(o,m)+pow(So,m))) + rho_m;
rho_v = rho;
E_v = (Eb-Em)*(pow(o,m)/(pow(o,m)+pow(So,m))) + Em;
nu_v = (nu_b-nu_m)*(pow(o,m)/(pow(o,m)+pow(So,m))) + nu_m;

mu = E_v/(2*(1+nu_v));
lambda = nu_v*E_v/((1+nu_v)*(1-2*nu_v));
muf= fvc::interpolate(mu);
lambdaf = fvc::interpolate(lambda); 





      int iCorr = 0;
      lduMatrix::solverPerformance solverPerf;
      scalar initialResidual = 1.0;
      scalar relativeResidual = 1.0;
      lduMatrix::debug = 0;

      
      
      do
      {
          DU.storePrevIter();

#         include "calculateDivDSigmaExp.H"
#         include "calculateDivDSigmaLargeStrainExp.H"
Info<< "Time2 = " << runTime.timeName() << nl << endl;
          //- Updated lagrangian momentum equation
          fvVectorMatrix DUEqn
              (
                  fvm::d2dt2(rho,DU)
                  ==
                  fvm::laplacian(2*muf + lambdaf, DU, "laplacian(DDU,DU)")
                  + divDSigmaExp
                  + divDSigmaLargeStrainExp
                  );



          if (solidInterfaceCorr)
          {
              solidInterfacePtr->correct(DUEqn);
          }

          solverPerf = DUEqn.solve();

          if (iCorr == 0)
          {
              initialResidual = solverPerf.initialResidual();
          }
	  
          DU.relax();
	  

          gradDU = fvc::grad(DU);
	  


#         include "calculateDEpsilonDSigma.H"
#         include "calculateRelativeResidual.H"



          Info << "\tTime " << runTime.value()
               << ", Corrector " << iCorr
               << ", Solving for " << DU.name()
               << " using " << solverPerf.solverName()
               << ", res = " << solverPerf.initialResidual()
               << ", rel res = " << relativeResidual
               << ", inner iters " << solverPerf.nIterations() << endl;
      }
      while
        (
            //solverPerf.initialResidual() > convergenceTolerance
            relativeResidual > convergenceTolerance
            && ++iCorr < nCorr
            );

      Info << nl << "Time " << runTime.value() << ", Solving for " << DU.name()
           << ", Initial residual = " << initialResidual
           << ", Final residual = " << solverPerf.initialResidual()
           << ", No outer iterations " << iCorr << endl;


#     include "moveMesh.H"

#     include "rotateFields.H"




#     include "writeFields.H"


      Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
          << "  ClockTime = " << runTime.elapsedClockTime() << " s"
          << endl;


epsilonHyd = tr(epsilon)/3;
epsilonEq = sqrt((2.0/3.0)*magSqr(dev(epsilon)));
epsilonSum = epsilonHyd+epsilonEq;

DEpsilonHyd = tr(DEpsilon)/3;
REpsilonHyd = DEpsilonHyd/(runTime.time().deltaT().value());
DEpsilonEq = sqrt((2.0/3.0)*magSqr(dev(DEpsilon)));
REpsilonEq = DEpsilonEq/(runTime.time().deltaT().value());
DEpsilonSum = 2*DEpsilonHyd+DEpsilonEq;

sigmaHyd = tr(sigma)/3;
sigmaEq = sqrt((3.0/2.0)*magSqr(dev(sigma)));
sigmaSum = sigmaHyd+sigmaEq;

DSigmaHyd = tr(DSigma)/3;
DSigmaEq = sqrt((3.0/2.0)*magSqr(dev(DSigma)));
DSigmaSum = 2*DSigmaHyd+DSigmaEq;

Info<< "Time = " << runTime.timeName() << nl << endl;



volTensorField F = I + gradDU;

volTensorField Finv = inv(F);

DEV = det(F)-1;
REV = DEV/(runTime.time().deltaT().value());
EV += DEV;

volVectorField v = DU/(runTime.time().deltaT()); 
surfaceScalarField phi = linearInterpolate(v) & mesh.Sf();

#       include "readSIMPLEControls.H"

	Da_v = (pow(So,m)/(pow(o,m)+pow(So,m)))*Da;
	Dh_v = (pow(So,m)/(pow(o,m)+pow(So,m)))*Dh;

          volTensorField    Da_v_T = Da*(Finv & Finv.T());

          volTensorField    Dh_v_T = Dh*(Finv & Finv.T());
	      


        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
             solve
            (
                fvm::ddt(a) 
		//-fvm::div(phi,a)
		//+(pow(So,m)/(pow(o,m)+pow(So,m)))*
		-a*fvc::div(phi)
		
	
	
-(pow(So,m)/(pow(o,m)+pow(So,m)))*((REpsilonHydOff+REV))*sigma_a

	
-(pow(So,m)/(pow(o,m)+pow(So,m)))*((REpsilonHydOff+REV))*sigma_o*o
		+(pow(So,m)/(pow(o,m)+pow(So,m)))*mu_a*a
		
		-(pow(So,m)/(pow(o,m)+pow(So,m)))*rho_a*a*a/h
		
		- fvm::laplacian(Da_v, a)
		
		
            );
	    
	    solve
            (
                fvm::ddt(h) 
		//-fvm::div(phi,h)
		-h*fvc::div(phi)
		-(pow(So,m)/(pow(o,m)+pow(So,m)))*sigma_h
		+(pow(So,m)/(pow(o,m)+pow(So,m)))*mu_h*h
		
		
		- (pow(So,m)/(pow(o,m)+pow(So,m)))*rho_h*a*a
		- fvm::laplacian(Dh_v, h)
		
		

            );
	



            solve
            (
                fvm::ddt(o)
		-o*fvc::div(phi)
               
-eta*(pow(So,m)/(pow(o,m)+pow(So,m)))*(pow((a*a/h),m)/(pow((a*a/h),m)+pow(aT,m)))*(pow(EV,m)/(pow(EV,m)+pow(leV,m)))
		
            ); 


	   
        }
	
aGrad = mag(fvc::grad(a));

oGrad = mag(fvc::grad(o));



UGrad = mag(fvc::grad(U));

	     
	     

        runTime.write();
	mesh.update();

DUGrad = mag(fvc::grad(DU));

    }

  Info<< "End\n" << endl;

  return(0);
}

// ************************************************************************* //
