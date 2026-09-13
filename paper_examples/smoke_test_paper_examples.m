% smoke_test_paper_examples.m — run from paper_examples/
setup_paths;
P = paper_optical_params();
fprintf('params ok nm=%.2f power=%.3f\n', P.nm, P.power);

shapes = {'dumbbell','fourBump','eightBump','neck','oblate75'};
for k = 1:numel(shapes)
    s = shapes{k};
    p = P.p_of(s);
    puse = min(p, 32);
    fprintf('building %s p=%d ... ', s, puse);
    particle = make_particle(s, puse, Point(0,0,0), [], [], [], P);
    fprintf('vol=%.4f area=%.4f\n', particle.cell.vol, particle.cell.surArea);
end

beam = CreateBeam(P.nm, P.f, P.NA, P.Ex0, P.Ey0, P.w0, P.Nphi, P.Nr, P.power);
particle = make_particle('dumbbell', 32, Point(0,0,0), [], [], [], P);
[F, T] = sum_force_torque(particle, beam, P.scat);
fprintf('F=[%.3e %.3e %.3e] N\n', F(1), F(2), F(3));
fprintf('T=[%.3e %.3e %.3e] N*m\n', T(1), T(2), T(3));

S = load(fullfile('..', 'otgo1.0.0', 'diffusion_tensors.mat'));
assert(isfield(S, 'D_dumbbell'), 'missing D_dumbbell');
D = S.D_dumbbell * 1e-7;
Lchol = chol(D, 'lower'); %#ok<NASGU>
fprintf('D chol ok size=%dx%d\n', size(D, 1), size(D, 2));

kB = PhysConst.kB;
Temp = 293;
dt = 2e-3;
Mlp = eye(3);
[Fl, Tl] = sum_force_torque(particle, beam, P.scat);
fp = Mlp * Fl;
tp = Mlp * Tl;
dq = ((D * dt) / (kB * Temp)) * [fp; tp] + sqrt(2 * dt) * (chol(D, 'lower') * randn(6, 1));
fprintf('one BD step |dq_trans|=%.3e\n', norm(dq(1:3)));

particle2 = make_particle('fourBump', 32, Point(0,0,0), [], [], [], P);
particle2 = particle2.xrotation(0.5);
[~, T2] = sum_force_torque(particle2, beam, P.scat);
fprintf('fourBump at 0.5 rad Tx=%.3e\n', T2(1));

% Displacement sample
c = Point(0.5e-6, 0, 0);
particle3 = make_particle('dumbbell', 32, c, particle.cell.or1, particle.cell.or2, particle.cell.or3, P);
[F3, ~] = sum_force_torque(particle3, beam, P.scat);
fprintf('dumbbell at x=0.5um Fx=%.3e N\n', F3(1));

fprintf('ALL SMOKE TESTS PASSED\n');
