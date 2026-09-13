function [F, T] = sum_force_torque(particle, beam, scat)
%SUM_FORCE_TORQUE Net lab-frame optical force [N] and torque [N·m].
%
%   [F,T] = sum_force_torque(particle, beam, scat)
%   F,T are 3x1 vectors [Fx;Fy;Fz] and [Tx;Ty;Tz].

r_vec = scattering(particle, beam, 0, scat - 1);
f = particle.force(r_vec, beam);
t = particle.torque(r_vec, beam);

F = [sum(f.Vx(isfinite(f.Vx))); ...
     sum(f.Vy(isfinite(f.Vy))); ...
     sum(f.Vz(isfinite(f.Vz)))];
T = [sum(t.Vx(isfinite(t.Vx))); ...
     sum(t.Vy(isfinite(t.Vy))); ...
     sum(t.Vz(isfinite(t.Vz)))];
end
