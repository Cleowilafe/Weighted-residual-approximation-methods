function u_at_x_05 = aprox_methods(N, l)
    % Função que implementa o método de Galerkin e outros métodos de resíduos ponderados
    %
    % Entradas:
    % N  - Número de termos da série
    % l  - Tipo de método de peso (1 = Galerkin, 2 = Least Squares, 3 = Método de Resíduo Ponderado)

    % Definições simbólicas
    syms x;
    a = sym('a', [1 N]); % Vetor de coeficientes simbólicos

    % Inicializando a função aproximadora
    u = 0;

    % Loop para adicionar os termos da série
    for k = 1:N
        u = u + a(k) * x^k * (1 - x); % Termo a_k * x^k * (1 - x)
    end

    % Definindo a função onde R é o resíduo
    f = diff(u, x, 2) - u + x; % Derivada segunda de u menos u mais x

    % Escolha do método de peso (w)
    w = escolher_metodo_w(u, f, a, l, N);

    % Resolvendo o sistema de equações
    eqs = int(f .* w, x, 0, 1) == 0;
    sol = solve(eqs, a);

    % Substituindo os valores de 'a_n' na função aproximadora
    for k = 1:N
        u = subs(u, a(k), sol(k)); % Substituindo a_k com o valor calculado de sol(k)
    end

    % Avaliando a função aproximadora em x = 0.5
    u_at_x_05 = double(subs(u, x, 0.5));

    % Exibindo os resultados
    disp('Valores de a_n (solução para a):');
    disp(sol);
    disp('Valor de u em x = 0.5:');
    disp(u_at_x_05);
end

function w = escolher_metodo_w(u, f, a, l, N)
    % Função para escolher o método de peso (w) baseado no valor de 'l'
    %
    % Entradas:
    % u  - Função aproximadora
    % f  - Resíduo (função)
    % a  - Parâmetro a
    % l  - Identificador do método de peso (1 = Galerkin, 2 = Least Squares, 3 = Resíduo Ponderado)
    % N  - Número de termos na série

    syms x;
    w = sym(zeros(1, N)); % Inicializando a função de peso

    switch l
        case 1  % Método de Galerkin
            for k = 1:N
                w(k) = diff(u, a(k)); % Derivada de u em relação a a_k
            end

        case 2  % Método dos Mínimos Quadrados
            for k = 1:N
                w(k) = diff(f, a(k)); % Derivada do resíduo em relação a a_k
            end

        case 3  % Método de Resíduo Ponderado
            for k = 1:N
                w(k) = x^k; % Função de peso genérica baseada em x^k
            end

        otherwise
            error('Método de peso desconhecido! Escolha entre 1 (Galerkin), 2 (Least Squares) ou 3 (Weighted Residual).');
    end
end

% Exemplo de uso:
N = 3; % Número de termos da série
l = 1; % Método de Galerkin (1 = Galerkin, 2 = Least Squares, 3 = Weighted Residual)

resultado = aprox_methods(N, l);
disp(['Resultado em x = 0.5: ', num2str(resultado)]);

