function [forces_torques, elapsedTime] = Solve(r,scat,power,par,data,N,bottom,top,order,varargin)

disp_scatterings = "yes";

if ~isempty(varargin)
    saveFilePath = varargin{1}; 
end

div = linspace(bottom,top,N);

switch order
    case 1
        
        tic;
        forces_torquesx = zeros(6,N);
        parfor i = 1:N
            particle = ParticleCSVCell(Point(div(i),par.cell.c.Y,par.cell.c.Z),par.cell.or1,par.cell.or2,par.cell.or3,"vol",-6,data,par.nm,par.np);
            r_vec = scattering(particle,r,0,scat-1,disp_scatterings);
            f = particle.force(r_vec,r);
            t = particle.torque(r_vec,r);
            
            force = [sum(f.Vx(isfinite(f.Vx)));
                         sum(f.Vy(isfinite(f.Vy)));
                         sum(f.Vz(isfinite(f.Vz)));
                        ];
            torque = [sum(t.Vx(isfinite(t.Vx)));
                          sum(t.Vy(isfinite(t.Vy)));
                          sum(t.Vz(isfinite(t.Vz)));
                         ];
        
            forces_torquesx(:,i) = [force(1); force(2); force(3); torque(1); torque(2); torque(3)];
        end
        
        elapsedTime = toc;
        fprintf('Total time X (N = %d, R = %d): %.2f seconds\n', N, numel(r), elapsedTime)

        if ~isempty(varargin)
            modifiedFilePath = sprintf('X(DISP)%d', N);
            modifiedFilePath = sprintf('%s(%dx%d)', modifiedFilePath, size(r,1),size(r,2));
            modifiedFilePath = sprintf('%s(scat%d)', modifiedFilePath,scat);
            modifiedFilePath = sprintf('%s(P%g)', modifiedFilePath, power);
            FilePath = fullfile(saveFilePath,modifiedFilePath);
            save(sprintf('%s%s',FilePath,'.mat'),'forces_torquesx');
        end

        forces_torques = forces_torquesx;

    case 2
        
        tic;
        forces_torquesy = zeros(6,N);
        parfor i = 1:N
            particle = ParticleCSVCell(Point(par.cell.c.X,div(i),par.cell.c.Z),par.cell.or1,par.cell.or2,par.cell.or3,"vol",-6,data,par.nm,par.np);
            r_vec = scattering(particle,r,0,scat-1,disp_scatterings);
            f = particle.force(r_vec,r);
            t = particle.torque(r_vec,r);
            
            force = [sum(f.Vx(isfinite(f.Vx)));
                         sum(f.Vy(isfinite(f.Vy)));
                         sum(f.Vz(isfinite(f.Vz)));
                        ];
            torque = [sum(t.Vx(isfinite(t.Vx)));
                          sum(t.Vy(isfinite(t.Vy)));
                          sum(t.Vz(isfinite(t.Vz)));
                         ];
        
            forces_torquesy(:,i) = [force(1); force(2); force(3); torque(1); torque(2); torque(3)];
        end
        
        elapsedTime = toc;
        fprintf('Total time Y (N = %d, R = %d): %.2f seconds\n', N, numel(r), elapsedTime)
        
        if ~isempty(varargin)
            modifiedFilePath = sprintf('Y(DISP)%d', N);
            modifiedFilePath = sprintf('%s(%dx%d)', modifiedFilePath, size(r,1),size(r,2));
            modifiedFilePath = sprintf('%s(scat%d)', modifiedFilePath,scat);
            modifiedFilePath = sprintf('%s(P%g)', modifiedFilePath, power);
            FilePath = fullfile(saveFilePath,modifiedFilePath);
            save(sprintf('%s%s',FilePath,'.mat'),'forces_torquesy');
        end
        
        forces_torques = forces_torquesy;

    case 3

        tic;
        forces_torquesz = zeros(6,N);
        parfor i = 1:N
            particle = ParticleCSVCell(Point(par.cell.c.X,par.cell.c.Y,div(i)),par.cell.or1,par.cell.or2,par.cell.or3,"vol",-6,data,par.nm,par.np);
            r_vec = scattering(particle,r,0,scat-1,disp_scatterings);
            f = particle.force(r_vec,r);
            t = particle.torque(r_vec,r);
            
            force = [sum(f.Vx(isfinite(f.Vx)));
                         sum(f.Vy(isfinite(f.Vy)));
                         sum(f.Vz(isfinite(f.Vz)));
                        ];
            torque = [sum(t.Vx(isfinite(t.Vx)));
                          sum(t.Vy(isfinite(t.Vy)));
                          sum(t.Vz(isfinite(t.Vz)));
                         ];
        
            forces_torquesz(:,i) = [force(1); force(2); force(3); torque(1); torque(2); torque(3)];
        end
        
        elapsedTime = toc;
        fprintf('Total time Z (N = %d, R = %d): %.2f seconds\n', N, numel(r), elapsedTime)
        
        if ~isempty(varargin)
            modifiedFilePath = sprintf('Z(DISP)%d', N);
            modifiedFilePath = sprintf('%s(%dx%d)', modifiedFilePath, size(r,1),size(r,2));
            modifiedFilePath = sprintf('%s(scat%d)', modifiedFilePath,scat);
            modifiedFilePath = sprintf('%s(P%g)', modifiedFilePath, power);
            FilePath = fullfile(saveFilePath,modifiedFilePath);
            save(sprintf('%s%s',FilePath,'.mat'),'forces_torquesz');
        end

        forces_torques = forces_torquesz;

    case 4

        tic;
        forces_torquesx_rot = zeros(6,N);
        parfor i = 1:N
            particle = ParticleCSVCell(Point(par.cell.c.X,par.cell.c.Y,par.cell.c.Z),par.cell.or1,par.cell.or2,par.cell.or3,"vol",-6,data,par.nm,par.np);
            particle = particle.xrotation(div(i))
        
            r_vec = scattering(particle,r,0,scat-1,disp_scatterings);
            f = particle.force(r_vec,r);
            t = particle.torque(r_vec,r);
            
            force = [sum(f.Vx(isfinite(f.Vx)));
                         sum(f.Vy(isfinite(f.Vy)));
                         sum(f.Vz(isfinite(f.Vz)));
                        ];
            torque = [sum(t.Vx(isfinite(t.Vx)));
                          sum(t.Vy(isfinite(t.Vy)));
                          sum(t.Vz(isfinite(t.Vz)));
                         ];
        
            forces_torquesx_rot(:,i) = [force(1); force(2); force(3); torque(1); torque(2); torque(3)];
        end
        
        elapsedTime = toc;
        fprintf('Total time X(ROT) (N = %d, R = %d): %.2f seconds\n', N, numel(r), elapsedTime)  

        if ~isempty(varargin)
            modifiedFilePath = sprintf('X(ROT)%d', N);
            modifiedFilePath = sprintf('%s(%dx%d)', modifiedFilePath, size(r,1),size(r,2));
            modifiedFilePath = sprintf('%s(scat%d)', modifiedFilePath,scat);
            modifiedFilePath = sprintf('%s(P%g)', modifiedFilePath, power);
            FilePath = fullfile(saveFilePath,modifiedFilePath);
            save(sprintf('%s%s',FilePath,'.mat'),'forces_torquesx_rot');
        end
        
        forces_torques = forces_torquesx_rot;

    case 5

        tic;
        forces_torquesy_rot = zeros(6,N);
        parfor i = 1:N
            particle = ParticleCSVCell(Point(par.cell.c.X,par.cell.c.Y,par.cell.c.Z),par.cell.or1,par.cell.or2,par.cell.or3,"vol",-6,data,par.nm,par.np);
            particle = particle.yrotation(div(i))
        
            r_vec = scattering(particle,r,0,scat-1,disp_scatterings);
            f = particle.force(r_vec,r);
            t = particle.torque(r_vec,r);
            
            force = [sum(f.Vx(isfinite(f.Vx)));
                         sum(f.Vy(isfinite(f.Vy)));
                         sum(f.Vz(isfinite(f.Vz)));
                        ];
            torque = [sum(t.Vx(isfinite(t.Vx)));
                          sum(t.Vy(isfinite(t.Vy)));
                          sum(t.Vz(isfinite(t.Vz)));
                         ];
        
            forces_torquesy_rot(:,i) = [force(1); force(2); force(3); torque(1); torque(2); torque(3)];
        end
        
        elapsedTime = toc;
        fprintf('Total time Y(ROT) (N = %d, R = %d): %.2f seconds\n', N, numel(r), elapsedTime)  

        if ~isempty(varargin)
            modifiedFilePath = sprintf('Y(ROT)%d', N);
            modifiedFilePath = sprintf('%s(%dx%d)', modifiedFilePath, size(r,1),size(r,2));
            modifiedFilePath = sprintf('%s(scat%d)', modifiedFilePath,scat);
            modifiedFilePath = sprintf('%s(P%g)', modifiedFilePath, power);
            FilePath = fullfile(saveFilePath,modifiedFilePath);
            save(sprintf('%s%s',FilePath,'.mat'),'forces_torquesy_rot');
        end

        forces_torques = forces_torquesy_rot;

    case 6

        tic;
        forces_torquesz_rot = zeros(6,N);
        parfor i = 1:N
            particle = ParticleCSVCell(Point(par.cell.c.X,par.cell.c.Y,par.cell.c.Z),par.cell.or1,par.cell.or2,par.cell.or3,"vol",-6,data,par.nm,par.np);
            particle = particle.zrotation(div(i))
        
            r_vec = scattering(particle,r,0,scat-1,disp_scatterings);
            f = particle.force(r_vec,r);
            t = particle.torque(r_vec,r);
            
            force = [sum(f.Vx(isfinite(f.Vx)));
                         sum(f.Vy(isfinite(f.Vy)));
                         sum(f.Vz(isfinite(f.Vz)));
                        ];
            torque = [sum(t.Vx(isfinite(t.Vx)));
                          sum(t.Vy(isfinite(t.Vy)));
                          sum(t.Vz(isfinite(t.Vz)));
                         ];
        
            forces_torquesz_rot(:,i) = [force(1); force(2); force(3); torque(1); torque(2); torque(3)];
        end
        
        elapsedTime = toc;
        fprintf('Total time Z(ROT) (N = %d, R = %d): %.2f seconds\n', N, numel(r), elapsedTime)  

        if ~isempty(varargin)
            modifiedFilePath = sprintf('Z(ROT)%d', N);
            modifiedFilePath = sprintf('%s(%dx%d)', modifiedFilePath, size(r,1),size(r,2));
            modifiedFilePath = sprintf('%s(scat%d)', modifiedFilePath,scat);
            modifiedFilePath = sprintf('%s(P%g)', modifiedFilePath, power);
            FilePath = fullfile(saveFilePath,modifiedFilePath);
            save(sprintf('%s%s',FilePath,'.mat'),'forces_torquesz_rot');
        end

        forces_torques = forces_torquesz_rot;
end

end

