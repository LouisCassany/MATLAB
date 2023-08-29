clc
close all
clearvars


% Répartition des points de fonctionement
discrete_time = [0 1 2:10:352 354:2:698 700:74:1440];

%% Plot des trajectoires d'états et des points de linéarisations

% Trajectoires des états du modèle pour un simple snénario de 50g sans
% bolus, données de référence issues de UVA/Padova
load sim_results.mat

for patient = 1
    
    % Vecteur temps et états issus d'UVA/Padova
    time50 = data.results(patient).time.signals.values;
    states50 = data.results(patient).state.signals.values(:,9:26);

    figure
    for i = 1:13
        plot(time50,states50(:,i));
        hold on
        plot(discrete_time,states50(discrete_time+1,i),'o','LineWidth',0.5)
        grid
        xlabel('Time (min)');
        ylabel(strcat('Model states'))
        legend("States trajectories",'Operatings points');
        title(strcat('Patient n°',num2str(patient)))
    end
    
end

%% Linéarisation autour de chacun des points de fonctionement

for patient = 1:11

    % Pour avoir les params_model et le basal du patient
    [~, params_model, ~, basal] = init_model(11+patient,1440, [360,15,50]);

    % Pour chaque point de fonctionement, linéarisation autour du point de
    % fonctionement
    for i=1:length(discrete_time)
        sys = get_ABCD(params_model, basal, states50(discrete_time(i)+1,:), states50(discrete_time(i)+1,:),i);
        A_arr(:,:,i) = sys.A;
        B_arr(:,:,i) = sys.B;
        C_arr(:,:,i) = sys.C;
        D_arr(:,:,i) = sys.D;
        disp(strcat('Patient n°',num2str(patient)))
        disp(strcat(num2str(i),'/',num2str(length(discrete_time))))
    end
    
    
    results(patient).A_arr = A_arr;
    results(patient).B_arr = B_arr;
    results(patient).C_arr = C_arr; 
    results(patient).D_arr = D_arr;
    
    patient
end


all_adults = results;

save('ABCD_all_adults','all_adults', 'discrete_time');