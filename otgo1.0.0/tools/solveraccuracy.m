function [acc, time, ftA, ftR] = solveraccuracy(r, scat, shp, c, p, nm, np, dis_dir, N, lims, varargin)

    if ~isempty(varargin)
        radius = varargin{1};
    else
        radius = 1*1e-6;
    end

    vol = 'vol';
    
    if strcmp(shp, 'sphere')
        %vol = (4/3)*pi*radius^3; % 10^-18 m^3
        shape_reference = ParticleSpherical(c,radius,nm,np);
        shape_actual = ParticleCELL(p, 'sphere', c, Vector(0,0,0,1,0,0), ...
            Vector(0,0,0,0,1,0), Vector(0,0,0,0,0,1), vol, -6, nm, np);
        
    elseif strcmp(shp, 'RBC')
        shape_reference = ParticleRBC(c,Vector(0,0,0,1,0,0),Vector(0,0,0,0,1,0),Vector(0,0,0,0,0,1),3.91*1e-6,0.81*1e-6,2.52*1e-6,2.76*1e-6,nm,np);
        shape_reference = shape_reference.xrotation(pi/2);
        shape_actual = ParticleCELL(p, 'RBC', c, Vector(0,0,0,1,0,0), ...
            Vector(0,0,0,0,1,0), Vector(0,0,0,0,0,1), vol, -6, nm, np);
        shape_actual = shape_actual.xrotation(pi/2);
    end

    or1 = shape_reference.rbc.or1;
    or2 = shape_reference.rbc.or2;
    or3 = shape_reference.rbc.or3;

    div = linspace(lims(1),lims(2),N);
    forces_torques_actual = zeros(6,N);
    forces_torques_reference = zeros(6,N);

    if dis_dir == 'X'
        tic;
        parfor i = 1:N
            particle_actual = ParticleCELL(p,shp,Point(div(i),c.Y,c.Z),shape_actual.cell.or1,shape_actual.cell.or2,shape_actual.cell.or3,vol,-6,nm,np);
            r_vec = scattering(particle_actual,r,0,scat-1);
            f = particle_actual.force(r_vec,r);
            t = particle_actual.torque(r_vec,r);
            
            force = [sum(f.Vx(isfinite(f.Vx)));
                         sum(f.Vy(isfinite(f.Vy)));
                         sum(f.Vz(isfinite(f.Vz)));
                        ];
            torque = [sum(t.Vx(isfinite(t.Vx)));
                          sum(t.Vy(isfinite(t.Vy)));
                          sum(t.Vz(isfinite(t.Vz)));
                         ];
        
            forces_torques_actual(:,i) = [force(1); force(2); force(3); torque(1); torque(2); torque(3)];
        end
        elapsedTime = toc;
        fprintf('Total time: %.2f seconds\n', elapsedTime)

        for i = 1:N
            if strcmp(shp, 'sphere')
                particle_reference = ParticleSpherical(Point(div(i),c.Y,c.Z),radius,nm,np);
            elseif strcmp(shp, 'RBC')
                particle_reference = ParticleRBC(Point(div(i),c.Y,c.Z),or1,or2,or3,3.91*1e-6,0.81*1e-6,2.52*1e-6,2.76*1e-6,nm,np);
            end
            r_vec = scattering(particle_reference,r,0,scat-1);
            f = particle_reference.force(r_vec,r);
            t = particle_reference.torque(r_vec,r);
            
            force = [sum(f.Vx(isfinite(f.Vx)));
                         sum(f.Vy(isfinite(f.Vy)));
                         sum(f.Vz(isfinite(f.Vz)));
                        ];
            torque = [sum(t.Vx(isfinite(t.Vx)));
                          sum(t.Vy(isfinite(t.Vy)));
                          sum(t.Vz(isfinite(t.Vz)));
                         ];
        
            forces_torques_reference(:,i) = [force(1); force(2); force(3); torque(1); torque(2); torque(3)];
        end

    elseif dis_dir == 'Y'
         tic;
         parfor i = 1:N
            particle_actual = ParticleCELL(p,shp,Point(c.X,div(i),c.Z),shape_actual.cell.or1,shape_actual.cell.or2,shape_actual.cell.or3,vol,-6,nm,np);
            r_vec = scattering(particle_actual,r,0,scat-1);
            f = particle_actual.force(r_vec,r);
            t = particle_actual.torque(r_vec,r);
            
            force = [sum(f.Vx(isfinite(f.Vx)));
                         sum(f.Vy(isfinite(f.Vy)));
                         sum(f.Vz(isfinite(f.Vz)));
                        ];
            torque = [sum(t.Vx(isfinite(t.Vx)));
                          sum(t.Vy(isfinite(t.Vy)));
                          sum(t.Vz(isfinite(t.Vz)));
                         ];
        
            forces_torques_actual(:,i) = [force(1); force(2); force(3); torque(1); torque(2); torque(3)];
         end
         elapsedTime = toc;
         fprintf('Total time: %.2f seconds\n', elapsedTime)

        for i = 1:N
            if strcmp(shp, 'sphere')
                particle_reference = ParticleSpherical(Point(c.X, div(i),c.Z),radius,nm,np);
            elseif strcmp(shp, 'RBC')
                particle_reference = ParticleRBC(Point(c.X, div(i),c.Z),or1,or2,or3,3.91*1e-6,0.81*1e-6,2.52*1e-6,2.76*1e-6,nm,np);
            end
            r_vec = scattering(particle_reference,r,0,scat-1);
            f = particle_reference.force(r_vec,r);
            t = particle_reference.torque(r_vec,r);
            
            force = [sum(f.Vx(isfinite(f.Vx)));
                         sum(f.Vy(isfinite(f.Vy)));
                         sum(f.Vz(isfinite(f.Vz)));
                        ];
            torque = [sum(t.Vx(isfinite(t.Vx)));
                          sum(t.Vy(isfinite(t.Vy)));
                          sum(t.Vz(isfinite(t.Vz)));
                         ];
        
            forces_torques_reference(:,i) = [force(1); force(2); force(3); torque(1); torque(2); torque(3)];
        end

    elseif dis_dir == 'Z'
        tic;
        parfor i = 1:N
            particle_actual = ParticleCELL(p,shp,Point(c.X,c.Y,div(i)),shape_actual.cell.or1,shape_actual.cell.or2,shape_actual.cell.or3,vol,-6,nm,np);
            r_vec = scattering(particle_actual,r,0,scat-1);
            f = particle_actual.force(r_vec,r);
            t = particle_actual.torque(r_vec,r);
            
            force = [sum(f.Vx(isfinite(f.Vx)));
                         sum(f.Vy(isfinite(f.Vy)));
                         sum(f.Vz(isfinite(f.Vz)));
                        ];
            torque = [sum(t.Vx(isfinite(t.Vx)));
                          sum(t.Vy(isfinite(t.Vy)));
                          sum(t.Vz(isfinite(t.Vz)));
                         ];
        
            forces_torques_actual(:,i) = [force(1); force(2); force(3); torque(1); torque(2); torque(3)];
        end
        elapsedTime = toc;
        fprintf('Total time: %.2f seconds\n', elapsedTime)

        for i = 1:N
            if strcmp(shp, 'sphere')
                particle_reference = ParticleSpherical(Point(c.X, c.Y, div(i)),radius,nm,np);
            elseif strcmp(shp, 'RBC')
                particle_reference = ParticleRBC(Point(c.X, c.Y,div(i)),or1,or2,or3,3.91*1e-6,0.81*1e-6,2.52*1e-6,2.76*1e-6,nm,np);
            end
            r_vec = scattering(particle_reference,r,0,scat-1);
            f = particle_reference.force(r_vec,r);
            t = particle_reference.torque(r_vec,r);
            
            force = [sum(f.Vx(isfinite(f.Vx)));
                         sum(f.Vy(isfinite(f.Vy)));
                         sum(f.Vz(isfinite(f.Vz)));
                        ];
            torque = [sum(t.Vx(isfinite(t.Vx)));
                          sum(t.Vy(isfinite(t.Vy)));
                          sum(t.Vz(isfinite(t.Vz)));
                         ];
        
            forces_torques_reference(:,i) = [force(1); force(2); force(3); torque(1); torque(2); torque(3)];
        end

    end

save('forces_torques_actual'); save('forces_torques_reference');
    
% forces
figure; hold on; box on;
plot(div,forces_torques_reference(1,:)*1e+12,'LineWidth',2)
plot(div,forces_torques_actual(1,:)*1e+12,'LineWidth',2)
xlabel(sprintf('%s coordinates (m)', dis_dir),'FontSize',18)
ylabel('X Force (pN)','FontSize',18)
legend('Reference X Force','Actual X Force','Location','best', 'Box', 'off')
set(gca,'linewidth',2)

figure; hold on; box on;
plot(div,forces_torques_reference(2,:)*1e+12,'LineWidth',2)
plot(div,forces_torques_actual(2,:)*1e+12,'LineWidth',2)
xlabel(sprintf('%s coordinates (m)', dis_dir),'FontSize',18)
ylabel('Y Force (pN)','FontSize',18)
legend('Reference Y Force','Actual Y Force','Location','best', 'Box', 'off')
set(gca,'linewidth',2)

figure; hold on; box on;
plot(div,forces_torques_reference(3,:)*1e+12,'LineWidth',2)
plot(div,forces_torques_actual(3,:)*1e+12,'LineWidth',2)
xlabel(sprintf('%s coordinates (m)', dis_dir),'FontSize',18)
ylabel('Z Force (pN)','FontSize',18)
legend('Reference Z Force','Actual Z Force','Location','best', 'Box', 'off')
set(gca,'linewidth',2)

% torques
figure; hold on; box on;
plot(div,forces_torques_reference(4,:),'LineWidth',2)
plot(div,forces_torques_actual(4,:),'LineWidth',2)
xlabel(sprintf('%s coordinates (m)', dis_dir),'FontSize',18)
ylabel('X Torque (Nm)','FontSize',18)
legend('Reference X Torque','Actual X Torque','Location','best', 'Box', 'off')
set(gca,'linewidth',2)

figure; hold on; box on;
plot(div,forces_torques_reference(5,:),'LineWidth',2)
plot(div,forces_torques_actual(5,:),'LineWidth',2)
xlabel(sprintf('%s coordinates (m)', dis_dir),'FontSize',18)
ylabel('Y Torque (Nm)','FontSize',18)
legend('Reference Y Torque','Actual Y Torque','Location','best', 'Box', 'off')
set(gca,'linewidth',2)

figure; hold on; box on;
plot(div,forces_torques_reference(6,:),'LineWidth',2)
plot(div,forces_torques_actual(6,:),'LineWidth',2)
xlabel(sprintf('%s coordinates (m)', dis_dir),'FontSize',18)
ylabel('Z Torque (Nm)','FontSize',18)
legend('Reference Z Torque','Actual Z Torque','Location','best', 'Box', 'off')
set(gca,'linewidth',2)

% differences
figure; hold on; box on;
plot(div,abs(forces_torques_reference(1,:)-forces_torques_actual(1,:))*1e+12,'LineWidth',2)
xlabel(sprintf('%s coordinates (m)', dis_dir),'FontSize',18)
ylabel('X Force (pN)','FontSize',18)
legend('Difference','Location','best', 'Box', 'off')
set(gca,'linewidth',2)

figure; hold on; box on;
plot(div,abs(forces_torques_reference(2,:)-forces_torques_actual(2,:))*1e+12,'LineWidth',2)
xlabel(sprintf('%s coordinates (m)', dis_dir),'FontSize',18)
ylabel('Y Force (pN)','FontSize',18)
legend('Difference','Location','best', 'Box', 'off')
set(gca,'linewidth',2)

figure; hold on; box on;
plot(div,abs(forces_torques_reference(3,:)-forces_torques_actual(3,:))*1e+12,'LineWidth',2)
xlabel(sprintf('%s coordinates (m)', dis_dir),'FontSize',18)
ylabel('Z Force (pN)','FontSize',18)
legend('Difference','Location','best', 'Box', 'off')
set(gca,'linewidth',2)

disp('*------------------------------------------------------------------------*')
fprintf('Mean Distance: %g\n', mean(mean(abs(forces_torques_reference-forces_torques_actual))));
fprintf('Total Distance: %g\n', sum(sum(abs(forces_torques_reference-forces_torques_actual))));
fprintf('Max Distance: %g\n', max(max(abs(forces_torques_reference-forces_torques_actual))));
fprintf('Root Mean Square Error: %g\n', sqrt(sum(sum(abs(forces_torques_reference-forces_torques_actual).^2)/numel(abs(forces_torques_reference-forces_torques_actual)))));

acc = [mean(mean(abs(forces_torques_reference-forces_torques_actual))), sum(sum(abs(forces_torques_reference-forces_torques_actual))), ...
       max(max(abs(forces_torques_reference-forces_torques_actual))), ...
       sqrt(sum(sum(abs(forces_torques_reference-forces_torques_actual).^2))/numel(abs(forces_torques_reference-forces_torques_actual)))];

time = elapsedTime;

ftA = forces_torques_actual;
ftR = forces_torques_reference;

end




