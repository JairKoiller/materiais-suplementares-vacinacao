
syms N S I R V real
syms mu beta phi sigma gamma theta psi positive

F = [mu*N-beta*S*I-mu*S-phi*S+theta*V+psi*R;...
     beta*S*I+sigma*beta*V*I-gamma*I-mu*I;...
     phi*S-mu*V-sigma*beta*V*I-theta*V;...
     gamma*I-mu*R-psi*R];
DF = jacobian(F, [S, I, R, V]);
%%
syms S I R V real

% para determinar os valores numéricos dos parametros, execute uma seção do
% script Paramters.m

Fnum = [mu*N-beta*S*I-mu*S-phi*S+theta*V+psi*R;...
     beta*S*I+sigma*beta*V*I-gamma*I-mu*I;...
     phi*S-mu*V-sigma*beta*V*I-theta*V;...
     gamma*I-mu*R-psi*R];

DFn = jacobian(Fnum, [S, I, V, R]);
sol  = solve(Fnum(1), Fnum(2), Fnum(3), Fnum(4), {S,I,V,R}, 'ReturnConditions', true);
for k = 1:length(sol.I)
    DFneq = double(subs(DFn, {S, I, V, R}, {sol.S(k), sol.I(k), sol.V(k), sol.R(k)}));
    [eigV, eigD] = eig(DFneq);
    disp(['** Equilibrium ', num2str(k), ' **'])
    disp(['     S=', num2str(double(sol.S(k))), ' I=', num2str(double(sol.I(k))), ' V=',...
        num2str(double(sol.V(k))), ' R=', num2str(double(sol.R(k)))])
    disp('     Eigenvalues')
    disp(['     ', num2str(transpose(diag(eigD)))])
end
%%
% Endemic equilibrium without vaccine
%I0 = 5.1703;
I0 = 0.1;
V0 = 0;
%S0 = 28.5714;
%R0 = 66.2582;
R0 = 0;

% some infections after reach the DFE
%I0 = 0.1;
%V0 = 90;
%S0 = 7.5;
%R0 = 0;
%yini = [parameters.N-I0-V0, I0, V0, R0];
yini = [parameters.N-(I0+V0+R0), I0, V0, R0];

[t, y] = ode45(@(t, y)(SIRV(y, parameters)),[0 3500],yini);

% --- Início da seção de plotagem melhorada ---
% Cria a figura com fundo branco e um tamanho fixo retangular (ideal para artigos)
figure('Color', 'w', 'Position', [100, 100, 800, 500]); 
hold on;

% Paleta de cores mais suave e profissional (RGB)
color_S = [0.4660 0.6740 0.1880]; % Verde
color_I = [0.8500 0.3250 0.0980]; % Laranja/Vermelho
color_V = [0.9290 0.6940 0.1250]; % Amarelo mostarda
color_R = [0.0000 0.4470 0.7410]; % Azul

% Espessura das linhas para melhor visualização
lw = 2.0;

% Plotagem (Comente as linhas das variáveis que não deseja mostrar)
%plot(t, y(:, 1), 'DisplayName', 'S (Suscetíveis)', 'Color', color_S, 'LineWidth', lw);
plot(t, y(:, 2), 'DisplayName', 'I (Infectados)', 'Color', color_I, 'LineWidth', lw);
plot(t, y(:, 3), 'DisplayName', 'V (Vacinados)', 'Color', color_V, 'LineWidth', lw);
plot(t, y(:, 4), 'DisplayName', 'R (Recuperados)', 'Color', color_R, 'LineWidth', lw);

% Formatação dos Eixos
xlabel('Tempo (dias)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Percentual de Indivíduos', 'FontSize', 12, 'FontWeight', 'bold');

% Estilização do Gráfico (Grade e Fontes)
set(gca, 'FontSize', 11, 'LineWidth', 1); % Ajusta tamanho dos números nos eixos
grid on;      % Grade principal
grid minor;   % Grade secundária (linhas pontilhadas)
box on;       % Fecha o contorno do gráfico

% Configuração da Legenda
legend('show', 'Location', 'best', 'FontSize', 11);
legend('boxoff'); % Remove a borda preta da legenda para um visual mais limpo

hold off;
% --- Fim da seção de plotagem melhorada ---
%disp(y(length(y),:))

%%
[i0m,v0m] = meshgrid(linspace(0.001,3,100), linspace(0,96,100));
Z = zeros(size(i0m));
cont = 0;
for lin = 1:size(Z,1)
    for col = 1:size(Z,2)
        cont = cont + 1; % Increment the counter for each iteration
        if mod(cont,100000) == 0
            disp(cont)
        end
        i0 = i0m(lin,col);
        v0 = v0m(lin,col);
        %r0 = 0;
        r0 = double(sol.R(1));
        yini = [parameters.N-i0-v0-r0, i0, v0, r0];
        [t, y] = ode45(@(t, y)(SIRV(y, parameters)),[0 3000],yini);
        If = y(length(y),2);
        %Z(lin,col) = If;
        if abs(If-2.3442) < 0.01
            Z(lin,col) = 3;
        elseif abs(If-0.1675) < 0.01
            Z(lin,col) = 2;
        elseif abs(If-0.0) < 0.01
            Z(lin,col) = 1;
        else
            Z(lin,col) = 4;
        end
        Z(lin,col) = If;
    end
end

%surf(i0m, v0m, Z)

%% --- 1. Configuração do Grid 3D com ZOOM ---
N = 100; 
res = 60; % Resolução (pode aumentar um pouco já que o volume é focado)

% --- INTERVALOS DE ZOOM ---
% I: Foco na região crítica perto de zero [0.001, 3]
vals_i = linspace(0, 3, res); 

% V: Alta cobertura vacinal [50, 100]
vals_v = linspace(50, 100, res);

% S: Matematicamente poderia ir até 100, mas como V >= 50, 
% o máximo S possível é 50. Ajustei para focar a resolução aqui.
vals_s = linspace(0, 50, res); 

[S_grid, I_grid, V_grid] = meshgrid(vals_s, vals_i, vals_v);

% Inicializa com NaN para ignorar pontos fora do simplexo no gráfico
Z = nan(size(S_grid)); 

%% --- 2. Simulação ---
total_points = numel(Z);
cont = 0;

fprintf('Iniciando simulação de %d pontos...\n', total_points);

for k = 1:size(Z,3)     % Loop V
    for j = 1:size(Z,2) % Loop S
        for i = 1:size(Z,1) % Loop I
            
            cont = cont + 1;
            
            s0 = S_grid(i,j,k);
            i0 = I_grid(i,j,k);
            v0 = V_grid(i,j,k);
            
            % --- RESTRIÇÃO DO SIMPLEXO ---
            if (s0 + i0 + v0) <= N
                r0 = N - s0 - i0 - v0;
                yini = [s0, i0, v0, r0];
                
                % Integração
                % Dica: Reduzi a tolerância (RelTol) para garantir precisão em I pequeno
                options = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);
                [~, y] = ode45(@(t, y) SIRV(y, parameters), [0 5000], yini, options);
                
                If = y(end, 2); 
                
                % --- CLASSIFICAÇÃO ---
                % Ajuste esses limiares conforme seus resultados observados
                if If < 0.01  % Tolerância para considerar "Zero"
                    Z(i,j,k) = 1; % Livre da Doença (Bacia Azul)
                elseif abs(If - 2.34) < 0.5 % Exemplo de Endêmico Alto
                    Z(i,j,k) = 3; % Endêmico (Bacia Vermelha)
                else
                    Z(i,j,k) = 2; % Região de transição ou outro equilíbrio
                end
                
            else
                % Fora do simplexo mantemos NaN
            end
        end
    end
end

fprintf('Simulação concluída.\n');

%% --- 3. Visualização Zoomada ---

figure('Color', 'w', 'Position', [100, 100, 1200, 600]);

% SUBPLOT 1: Scatter das Bacias
subplot(1, 2, 1); hold on;
title('Zoom nas Bacias de Atração');

% Pontos da Bacia Livre da Doença (Azul)
idx1 = Z == 1;
scatter3(S_grid(idx1), I_grid(idx1), V_grid(idx1), 20, 'b', 'filled', ...
    'MarkerFaceAlpha', 0.4, 'DisplayName', 'Livre da Doença');

% Pontos da Bacia Endêmica (Vermelho)
idx3 = Z == 3;
scatter3(S_grid(idx3), I_grid(idx3), V_grid(idx3), 20, 'r', 'filled', ...
    'MarkerFaceAlpha', 0.4, 'DisplayName', 'Endêmico');

xlabel('Suscetíveis'); 
ylabel('Infectados (Zoom)'); 
zlabel('Vacinados (Zoom)');
view(3); grid on; 
legend('Location', 'best');

% *** TRUQUE DE VISUALIZAÇÃO ***
% Como I vai de 0 a 3 e V de 50 a 100, o gráfico fica achatado.
% pbaspect força a caixa do gráfico a ser um cubo visualmente
pbaspect([1 1 1]); 


% SUBPLOT 2: Superfície Separatriz
subplot(1, 2, 2); hold on;
title('Superfície Separatriz (Fronteira)');

% Para isosurface funcionar bem com NaNs, as vezes precisamos de cuidado
% Se Z tiver NaNs, o isosurface pode ficar "buraquento".
% Aqui tentamos plotar direto.
val_fronteira = 2; % Valor entre as classes 1 e 3

p = patch(isosurface(S_grid, I_grid, V_grid, Z, val_fronteira));
p.FaceColor = [0.5 0.5 0.5];
p.EdgeColor = 'none';
p.FaceAlpha = 0.7;

isonormals(S_grid, I_grid, V_grid, Z, p);
camlight; lighting gouraud;

xlabel('Suscetíveis'); ylabel('Infectados'); zlabel('Vacinados');
grid on; view(120, 30);
pbaspect([1 1 1]); % Mantém a proporção visual cúbica

% Recriar os grids planos baseados nos vetores de zoom
[I_plane, V_plane] = meshgrid(vals_i, vals_v);

% Calcular S baseado na restrição do simplexo
S_plot = N - I_plane - V_plane;

% Criar uma máscara para o que é inválido
% Inválido se S < 0 (fora do simplexo físico) 
% OU se S > max(vals_s) (fora do nosso cubo de zoom visual)
mask_invalida = (S_plot < 0) | (S_plot > max(vals_s));

% Aplicar a máscara: transformamos os pontos ruins em NaN
S_plot(mask_invalida) = NaN;
% Se quiser garantir, pode aplicar a mesma máscara em I e V, mas geralmente S basta
% I_plane(mask_invalida) = NaN; 
% V_plane(mask_invalida) = NaN;

% Agora o surf recebe matrizes completas (com buracos de NaN)
surf(S_plot, I_plane, V_plane, ...
    'FaceColor', 'g', ...
    'FaceAlpha', 0.1, ...
    'EdgeColor', 'none', ...
    'DisplayName', 'Plano S+I+V=N');

% fig2plotly(gcf, ...
%     'filename', 'Bacias_Atracao_3D', ...
%     'offline', true, ...
%     'open', false, ...
%     'strip', false); 

%fprintf('Arquivo HTML gerado! Procure por "Bacias_Atracao_3D.html" no painel à esquerda.\n');