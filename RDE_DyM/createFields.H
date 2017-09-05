 

    Info<< "Reading incremental displacement field DU\n" << endl;
    volVectorField DU
    (
        IOobject
        (
            "DU",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volTensorField gradDU = fvc::grad(DU);

    Info<< "Reading accumulated displacement field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );
    
/*    surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(DU) & mesh.Sf()
);*/

volScalarField DUGrad
(IOobject
("DUGrad",
runTime.timeName(),
mesh,
IOobject::READ_IF_PRESENT,
IOobject::AUTO_WRITE
),
mag(fvc::grad(DU))
);



volScalarField UGrad
(IOobject
("UGrad",
runTime.timeName(),
mesh,
IOobject::READ_IF_PRESENT,
IOobject::AUTO_WRITE
),
mag(fvc::grad(U))
);


    Info << "Reading accumulated strain field epsilon\n" << endl;
    volSymmTensorField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );



    volSymmTensorField DEpsilon
    (
        IOobject
        (
            "DEpsilon",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

    Info << "Reading accumulated stress field sigma\n" << endl;
    volSymmTensorField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    );


    Info << "Reading incremental stress field DSigma\n" << endl;
    volSymmTensorField DSigma
    (
        IOobject
        (
            "DSigma",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    );

    //- explicit terms in the momentum equation
    volVectorField divDSigmaExp
	(
	 IOobject
	 (
	  "divDSigmaExp",
	  runTime.timeName(),
	  mesh,
	  IOobject::NO_READ,
	  IOobject::NO_WRITE
	  ),
	 mesh,
 	 dimensionedVector("zero", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector::zero)
	 );

     volVectorField divDSigmaLargeStrainExp
	(
	 IOobject
	 (
	  "divDSigmaLargeStrainExp",
	  runTime.timeName(),
	  mesh,
	  IOobject::NO_READ,
	  IOobject::NO_WRITE
	  ),
	 mesh,
 	 dimensionedVector("zero", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector::zero)
	 );

    constitutiveModel rheology(sigma, DU);

    volScalarField rho = rheology.rho();

 
volScalarField rho_v
    (
        IOobject
        (
            "rho_v",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


	
    volScalarField mu = rheology.mu();
    volScalarField lambda = rheology.lambda();
    surfaceScalarField muf = fvc::interpolate(rheology.mu());
    surfaceScalarField lambdaf = fvc::interpolate(rheology.lambda());

    surfaceVectorField n = mesh.Sf()/mesh.magSf();





    volScalarField E_v
    (
        IOobject
        (
            "E_v",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    volScalarField nu_v
    (
        IOobject
        (
            "nu_v",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    









    Info<< "Reading field a\n" << endl;

    volScalarField a
    (
        IOobject
        (
            "a",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
volScalarField aGrad
(IOobject
("aGrad",
runTime.timeName(),
mesh,
IOobject::READ_IF_PRESENT,
IOobject::AUTO_WRITE
),
mag(fvc::grad(a))
);



  Info<< "Reading field h\n" << endl;

    volScalarField h
    (
        IOobject
        (
            "h",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


Info<< "Reading field o\n" << endl;

    volScalarField o
    (
        IOobject
        (
            "o",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
volScalarField oGrad
(IOobject
("oGrad",
runTime.timeName(),
mesh,
IOobject::READ_IF_PRESENT,
IOobject::AUTO_WRITE
),
mag(fvc::grad(o))
);


 volScalarField Da_v
    (
        IOobject
        (
            "Da_v",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    
    );
    

volScalarField Dh_v
    (
        IOobject
        (
            "Dh_v",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    
    );
    
    



    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

Info<< "Reading sigma_a\n" << endl;

    dimensionedScalar sigma_a
    (
        transportProperties.lookup("sigma_a")
    );

Info<< "Reading sigma_h\n" << endl;

    dimensionedScalar sigma_h
    (
        transportProperties.lookup("sigma_h")
    );

Info<< "Reading sigma_o\n" << endl;

    dimensionedScalar sigma_o
    (
        transportProperties.lookup("sigma_o")
    );


Info<< "Reading mu_a\n" << endl;

    dimensionedScalar mu_a
    (
        transportProperties.lookup("mu_a")
    );

Info<< "Reading mu_h\n" << endl;

    dimensionedScalar mu_h
    (
        transportProperties.lookup("mu_h")
    );




Info<< "Reading rho_a\n" << endl;

    dimensionedScalar rho_a
    (
        transportProperties.lookup("rho_a")
    );

Info<< "Reading rho_h\n" << endl;

    dimensionedScalar rho_h
    (
        transportProperties.lookup("rho_h")
    );



 Info<< "Reading Da\n" << endl;

    dimensionedScalar Da
    (
        transportProperties.lookup("Da")
    );

 Info<< "Reading Dh\n" << endl;

    dimensionedScalar Dh
    (
        transportProperties.lookup("Dh")
    );



 Info<< "Reading aT\n" << endl;

    dimensionedScalar aT
    (
        transportProperties.lookup("aT")
    );
 
 Info<< "Reading hl\n" << endl;

    dimensionedScalar hl
    (
        transportProperties.lookup("hl")
    );

   
    Info<< "Reading So\n" << endl;

    dimensionedScalar So
    (
        transportProperties.lookup("So")
    );

 Info<< "Reading m\n" << endl;

    dimensionedScalar m
    (
        transportProperties.lookup("m")
    );

Info<< "Reading r\n" << endl;

    dimensionedScalar r
    (
        transportProperties.lookup("r")
    );


Info<< "Reading eta\n" << endl;

    dimensionedScalar eta
    (
        transportProperties.lookup("eta")
    );



    Info<< "Reading Eb\n" << endl;

    dimensionedScalar Eb
    (
        transportProperties.lookup("Eb")
    );

    Info<< "Reading Em\n" << endl;

    dimensionedScalar Em
    (
        transportProperties.lookup("Em")
    );

    Info<< "Reading nu_b\n" << endl;

    dimensionedScalar nu_b
    (
        transportProperties.lookup("nu_b")
    );

    Info<< "Reading nu_m\n" << endl;

    dimensionedScalar nu_m
    (
        transportProperties.lookup("nu_m")
    );


    Info<< "Reading rho_b\n" << endl;

    dimensionedScalar rho_b
    (
        transportProperties.lookup("rho_b")
    );


    Info<< "Reading rho_m\n" << endl;

    dimensionedScalar rho_m
    (
        transportProperties.lookup("rho_m")
    );


    Info<< "Reading leV\n" << endl;

    dimensionedScalar leV
    (
        transportProperties.lookup("leV")
    );

    Info<< "Reading leS\n" << endl;

    dimensionedScalar leS
    (
        transportProperties.lookup("leS")
    );


    Info<< "Reading DEpsilonHydOff\n" << endl;

    dimensionedScalar DEpsilonHydOff
    (
        transportProperties.lookup("DEpsilonHydOff")
    );


    Info<< "Reading REpsilonHydOff\n" << endl;

    dimensionedScalar REpsilonHydOff
    (
        transportProperties.lookup("REpsilonHydOff")
    );


volScalarField epsilonHyd
      (
       IOobject
       (
	"epsilonHyd",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
       tr(epsilon)/3
       );

volScalarField epsilonEq
      (
       IOobject
       (
	"epsilonEq",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
       sqrt((2.0/3.0)*magSqr(dev(epsilon)))
       );


volScalarField epsilonSum
      (
       IOobject
       (
	"epsilonSum",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
       epsilonEq+epsilonHyd
       );

volScalarField DEpsilonHyd
      (
       IOobject
       (
        "DEpsilonHyd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
       tr(DEpsilon)/3
       );
       
       volScalarField REpsilonHyd
      (
       IOobject
       (
        "REpsilonHyd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
       DEpsilonHyd/(runTime.time().deltaT().value())
       );

volScalarField DEpsilonEq
      (
       IOobject
       (
        "DEpsilonEq",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
       sqrt((2.0/3.0)*magSqr(dev(DEpsilon)))
       );

volScalarField REpsilonEq
      (
       IOobject
       (
        "REpsilonEq",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
       DEpsilonEq/(runTime.time().deltaT().value())
       );

volScalarField DEpsilonSum
      (
       IOobject
       (
        "DEpsilonSum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
       DEpsilonEq+2*DEpsilonHyd
       );


volScalarField sigmaHyd
      (
       IOobject
       (
	"sigmaHyd",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
       tr(sigma)/3
       );

volScalarField sigmaEq
      (
       IOobject
       (
	"sigmaEq",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
       sqrt((3.0/2.0)*magSqr(dev(sigma)))
       );


volScalarField sigmaSum
      (
       IOobject
       (
	"sigmaSum",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
	),
       sigmaEq+sigmaHyd
       );

volScalarField DSigmaHyd
      (
       IOobject
       (
        "DSigmaHyd",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
       tr(DSigma)/3
       );

volScalarField DSigmaEq
      (
       IOobject
       (
        "DSigmaEq",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
       sqrt((3.0/2.0)*magSqr(dev(DSigma)))
       );


volScalarField DSigmaSum
      (
       IOobject
       (
        "DSigmaSum",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
       DSigmaEq+2*DSigmaHyd
       );

       volScalarField cellVolume
      (
       IOobject
       (
        "cellVolume",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
        ),
       mesh,
       dimensionedScalar("zero",dimVolume,0.0)
       );
       cellVolume.internalField() = mesh.V();
       
       volScalarField dV
      (
       IOobject
       (
        "dV",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
        ),
       mesh,
       dimensionedScalar("zero",dimensionSet(0, 3, -1, 0, 0, 0, 0),0.0)
       );
       
        volScalarField DEV
      (
       IOobject
       (
        "DEV",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
        ),
       mesh,
       dimensionedScalar("zero",dimensionSet(0, 0, 0, 0, 0, 0, 0),0.0)
              );
	      
	volScalarField REV
      (
       IOobject
       (
        "REV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
        ),
       DEV/(runTime.time().deltaT().value())
       );	      
       
        volScalarField EV
      (
       IOobject
       (
        "EV",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
        ),
       mesh,
       dimensionedScalar("zero",dimensionSet(0, 0, 0, 0, 0, 0, 0),0.0)
              );
