function BlochAllObs(povms)
    nr_parties = size(povms,2);
    nr_inputs_per_party = [];
    nr_outputs_per_party = [];
    for party=1:nr_parties
        inps = size(povms{party},2);
        outs = size(povms{party}{1},1);
        nr_inputs_per_party = [nr_inputs_per_party, inps];
    end
    
    bloch_components = {{}};
    string_names = {{}};
    for party=1:nr_parties
       party_char = char(party+'A'-1);
       for x=1:nr_inputs_per_party(party)
           party_in = string(x);
           string_names{party}{x} = strcat(party_char, party_in);
           allbloch = BlochComponents(povms{party}{x}{1}-povms{party}{x}{2});
           bloch_components{party}{x} = allbloch(2:end); % remove id
       end
    end
    
    
    figure(1)
    hold on
    
    [Xs, Yx, Zx] = sphere(25);
    mySphere = surf( Xs, Yx, Zx );
    axis equal
    shading interp
    xlabel('x')
    ylabel('y')
    zlabel('z')
    mySphere.FaceAlpha = 0.25;
    
    line([-1 1], [0 0], [0 0])
    line([0 0], [-1 1], [0 0])
    line([0 0], [0 0], [-1 1])
    
    text(0, 0, 1.1, '$\left| 0 \right>$', 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'Center')
    text(1.1, 0, 0, '$\left| + \right>$', 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'Center')
    text(-1.1, 0, 0, '$\left| - \right>$', 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'Center')    
    text(0, 0, -1.1, '$\left| 1 \right>$', 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'Center')
    

    view([60 12])
    party_colors = {'r','b','g'};
    past_vectors = [[0;0;0]];
    for party=1:nr_parties
        for x=1:nr_inputs_per_party(party)
            p1 = [0 0 0];
            p2 = bloch_components{party}{x} ;
            dp = p2-p1;
%             past_vectors = [past_vectors, p2'];
%             flag = false;
%             while flag == false
%                 p2 = p2 + [0.1,0.1,0.1];
%                 for vec=1:size(past_vectors,2)
%                     if past_vectors(:,vec) == p2
%                     flag = false;
%                     break;
%                     end
%                 end
%                 disp(past_vectors);
%             end
            text(p2(1)+0.1, p2(2)+0.1, p2(3)+0.1, string_names{party}{x}, 'FontSize', 14, 'HorizontalAlignment', 'Center')
            
            quiver3(p1(1),p1(2),p1(3),dp(1),dp(2),dp(3),'LineWidth',3,'DisplayName',string_names{party}{x},'Color',party_colors{party})

        end
    end
    
    legend()
    
end

