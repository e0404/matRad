function test_suite = test_rayTracer

    test_functions=localfunctions();

    initTestSuite;
    
    function test_siddeonRayTracer

    % test funcion with dummy nummerical example

    cubes{1} = ones([2,2,2]);
    cubes{2} = cubes{1};
    cubes{2}(:,:,2) = [2,2; 2,2];

    resolution.x = 1;
    resolution.y = 1;
    resolution.z = 1;

    isocenter   = [0,0,0];
    sourcePoint = [1.5, 1.5, -4];
    targetPoint = [ 2.5,  2.5, 6];

   
    [alphas,l,rho,d12,ix] = matRad_siddonRayTracer(isocenter,...
            resolution, ...
            sourcePoint, ...
            targetPoint, ...
            cubes);

    % test Output types
    assertTrue(isvector(alphas));
    assertTrue(isvector(l));
    assertTrue(iscell(rho));
    assertTrue(isfloat(d12));  
    assertTrue(isvector(ix));

    % test numerical Output
    entryPoints = [1.95,1.95,0.5;
                    2.05,2.05,1.5;
                    2.15,2.15,2.5];
    entryPoints = entryPoints - sourcePoint;
    alphasNum   = sqrt(sum(entryPoints.^2,2))./sqrt(102);
    lNum        = [sqrt(102)/10,sqrt(102)/10];
    rhoNum{1}   = [1,1];
    rhoNum{2}   = [1,2];
    d12Num      = sqrt(102);
    ixNum       = [4,8];

    assertElementsAlmostEqual(alphasNum',alphas)
    assertElementsAlmostEqual(lNum,l)
    assertElementsAlmostEqual(rhoNum{1},rho{1})
    assertElementsAlmostEqual(rhoNum{2},rho{2})
    assertElementsAlmostEqual(d12Num,d12)
    assertEqual(ixNum,ix)

    function test_2DCube

        cubes{1} = ones([2,2]);
        cubes{2} = cubes{1};
        cubes{2}(:,2) = [2,2];
    
        resolution.x = 1;
        resolution.y = 1;
        resolution.z = 1;
    
        isocenter   = [0,0,0];
        sourcePoint = [1.5, 1.5, -4];
        targetPoint = [ 2.5,  2.5, 6];
    
       
        [alphas,l,rho,d12,ix] = matRad_siddonRayTracer(isocenter,...
                resolution, ...
                sourcePoint, ...
                targetPoint, ...
                cubes);

        % test Output types
        assertTrue(isvector(alphas));
        assertTrue(isvector(l));
        assertTrue(iscell(rho));
        assertTrue(isfloat(d12));  
        assertTrue(isvector(ix));

        % test numerical Output
        entryPoints = [1.95,1.95,0.5;
                        2.05,2.05,1.5;];
        entryPoints = entryPoints - sourcePoint;
        alphasNum   = sqrt(sum(entryPoints.^2,2))./sqrt(102);
        lNum        = [sqrt(102)/10];
        rhoNum{1}   = [1];
        rhoNum{2}   = [2];
        d12Num      = sqrt(102);
        ixNum       = [4];

        assertElementsAlmostEqual(alphasNum',alphas)
        assertElementsAlmostEqual(lNum,l)
        assertElementsAlmostEqual(rhoNum{1},rho{1})
        assertElementsAlmostEqual(rhoNum{2},rho{2})
        assertElementsAlmostEqual(d12Num,d12)
        assertEqual(ixNum,ix)

    function test_rayDoesNotHitCT
        
        cubes{1} = ones([2,2,2]);

        resolution.x = 1;
        resolution.y = 1;
        resolution.z = 1;

        isocenter   = [0,0,0];
        sourcePoint = [1.5, 1.5 -4];
        targetPoint = [ 10, 10, 6];

    
        [alphas,l,rho,d12,ix] = matRad_siddonRayTracer(isocenter,...
                resolution, ...
                sourcePoint, ...
                targetPoint, ...
                cubes);


        % test numerical Output
        alphasNum   = [];
        lNum        = [];
        rhoNum{1}   = [];
        d12Num      = norm(sourcePoint - targetPoint);
        ixNum       = [];

        
        assertElementsAlmostEqual(alphasNum',alphas)
        assertElementsAlmostEqual(lNum,l)
        assertElementsAlmostEqual(rhoNum{1},rho{1})
        assertElementsAlmostEqual(d12Num,d12)
        assertEqual(ixNum,ix)

        sourcePoint = [10, 10 -4];
        targetPoint = [ 10, 10, 6];
        d12Num      = 10;

    
        [alphas,l,rho,d12,ix] = matRad_siddonRayTracer(isocenter,...
                resolution, ...
                sourcePoint, ...
                targetPoint, ...
                cubes);
        
        assertElementsAlmostEqual(alphasNum',alphas)
        assertElementsAlmostEqual(lNum,l)
        assertElementsAlmostEqual(rhoNum{1},rho{1})
        assertElementsAlmostEqual(d12Num,d12)
        assertEqual(ixNum,ix)

    function test_rayHitsAtBoundary

        cubes{1} = ones([2,2,2]);
        cubes{2} = cubes{1};
        cubes{2}(:,:,2) = [2,2; 2,2];
    
        resolution.x = 1;
        resolution.y = 1;
        resolution.z = 1;
    
        isocenter   = [0,0,0];
        sourcePoint = [2.5, 2.5, -4];
        targetPoint = [2.5, 2.5, 6];
    
       
        [alphas,l,rho,d12,ix] = matRad_siddonRayTracer(isocenter,...
                resolution, ...
                sourcePoint, ...
                targetPoint, ...
                cubes);
        
        % test Output types
        assertTrue(isvector(alphas));
        assertTrue(isvector(l));
        assertTrue(iscell(rho));
        assertTrue(isfloat(d12));  
        assertTrue(isvector(ix));

        % test numerical Output
        entryPoints = [2.5, 2.5, 0.5;
                       2.5, 2.5, 1.5;
                       2.5, 2.5, 2.5];
        entryPoints = entryPoints - sourcePoint;
        alphasNum   = sqrt(sum(entryPoints.^2,2))./10;
        lNum        = [1,1];
        rhoNum{1}   = [1,1];
        rhoNum{2}   = [1,2];
        d12Num      = 10;
        ixNum       = [4,8];

        assertElementsAlmostEqual(alphasNum',alphas)
        assertElementsAlmostEqual(lNum,l)
        assertElementsAlmostEqual(rhoNum{1},rho{1})
        assertElementsAlmostEqual(rhoNum{2},rho{2})
        assertElementsAlmostEqual(d12Num,d12)
        assertEqual(ixNum,ix)

    
        function test_rayHitsAtCorner

            cubes{1} = ones([2,2,2]);
            cubes{2} = cubes{1};
            cubes{2}(:,:,2) = [2,2; 2,2];
        
            resolution.x = 1;
            resolution.y = 1;
            resolution.z = 1;
        
            isocenter   = [0,0,0];
            sourcePoint = [1.5, 1.5, -4];
            targetPoint = [3.5, 3.5, 5];
        
           
            [alphas,l,rho,d12,ix] = matRad_siddonRayTracer(isocenter,...
                    resolution, ...
                    sourcePoint, ...
                    targetPoint, ...
                    cubes);
        
    
            % test numerical Output
            alphasNum   = 0.5;
            ixNum       = 1:0;
            lNum        = [];
            rhoNum{1}   = cubes{1}(ixNum);
            rhoNum{2}   = cubes{2}(ixNum);
            d12Num      = norm(sourcePoint - targetPoint);

    
            assertElementsAlmostEqual(alphasNum',alphas)
            assertElementsAlmostEqual(lNum,l)
            assertElementsAlmostEqual(rhoNum{1},rho{1})
            assertElementsAlmostEqual(rhoNum{2},rho{2})
            assertElementsAlmostEqual(d12Num,d12)
            assertEqual(ixNum,ix)

        
