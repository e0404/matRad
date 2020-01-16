function mask = slabGeometry(vars, cubeDim)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% takes the shape and the size of the geometry and produces the
% corresponding mask cube for it.
%
%   call:
%         mask = slabGeometry(vars, cubeDim)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geo = vars.geoSize;
slab_loc = vars.slab_loc;

mask = zeros(cubeDim);

switch vars.geoShape
    
    case 'Rectangle'
        
        for i = -geo(1):geo(1)
            for j = -geo(2) : geo(2)
                for k = -geo(3):geo(3)
                    ix = slab_loc + [i, j, k];
                    if(nnz(ix>0) == 3)
                        mask(ix(1),ix(2),ix(3)) = 1;
                    end
                end
            end
        end
        
    case 'Ellipsoid'
        
        for i = -geo(1):geo(1)
            for j = -geo(2):geo(2)
                for k = -geo(3):geo(3)
                    ix = slab_loc + [i,j,k];
                    if (sum(([i,j,k].^2) ./ (geo.^2)) <= 1 && nnz(ix>0) == 3)
                        mask(ix(1),ix(2),ix(3)) = 1;
                    end                       
                end
            end
        end
    case '2DPyramid'
        
        for i = -geo(1):geo(1)
            for j = -geo(2):geo(2)
                for k = -geo(3):geo(3)
                    ix = slab_loc + [i,j,k];
                    if geo(2)
                        coef = (1 - sign(j)*j/geo(2));
                        assert(coef <= 1);
                    else
                        coef = 1;
                    end
                    
                    if (sign(i)*i <= coef * geo(1)  && nnz(ix>0) == 3 )
                        mask(ix(1),ix(2),ix(3)) = 1;
                    end                       
                end
            end
        end
        
    case 'Pyramid'
        
        for i = -geo(1):geo(1)
            for j = -geo(2):geo(2)
                for k = -geo(3):geo(3)
                    ix = slab_loc + [i,j,k];
                    if geo(2)
                        coef = (1 - sign(j)*j/geo(2));
                        assert(coef <= 1);
                    else
                        coef = 1;
                    end
                    
                    if (sign(i)*i <= coef * geo(1) && sign(k)*k <= coef * geo(3) && nnz(ix>0) == 3)
                        mask(ix(1),ix(2),ix(3)) = 1;
                    end                       
                end
            end
        end
        
    otherwise
        
            error('the slab geometry is not defined')
            
end

% 
%     center = slab_loc - [0 round(geo(2)/2) 0];
%                     currPos = center + [i j z];
                             
end