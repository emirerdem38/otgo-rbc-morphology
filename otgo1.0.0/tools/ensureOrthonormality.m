function Mlp = ensureOrthonormality(Mlp)
    % Perform Gram-Schmidt orthonormalization
    Mlp(:, 1) = Mlp(:, 1) / norm(Mlp(:, 1));
    Mlp(:, 2) = Mlp(:, 2) - dot(Mlp(:, 2), Mlp(:, 1)) * Mlp(:, 1);
    Mlp(:, 2) = Mlp(:, 2) / norm(Mlp(:, 2));
    Mlp(:, 3) = Mlp(:, 3) - dot(Mlp(:, 3), Mlp(:, 1)) * Mlp(:, 1) - dot(Mlp(:, 3), Mlp(:, 2)) * Mlp(:, 2);
    Mlp(:, 3) = Mlp(:, 3) / norm(Mlp(:, 3));
end