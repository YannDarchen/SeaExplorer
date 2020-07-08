dir27 = 'C:\Users\yannd\Documents\';

nav_all_file = {'sea027.95.gli.sub'};

    
% SEA27 nav file

for i = 1

   datanav = importdata([dir27 nav_all_file{1} '.all'],';'); 

   index_txt = max(find(~cellfun(@isempty,datanav.textdata(2,:))));

   header_nav = datanav.textdata(1,:);

    % Extract data

    i = 1;

    for i = 1:length(header_nav)

        if i > index_txt % .data

            index = i - index_txt;

            eval([strcat('SEA27_NAV_',header_nav{i}) '= datanav.data(:,index);'])

        else % .textdata

            index = i;

            if i == 3 % DATE

                eval([strcat('SEA27_NAV_',header_nav{i}) '= char(datanav.textdata((2:end),index));'])

            else

                eval([strcat('SEA27_NAV_',header_nav{i}) '= str2num(char(datanav.textdata((2:end),index)));'])

            end

        end 

    end

   

    

for i = 1:length(SEA27_NAV_Timestamp)

    SEA27_NAV_TIME(i) = datenum(SEA27_NAV_Timestamp(i,:),'dd/mm/yyyy HH:MM:SS');

end


end