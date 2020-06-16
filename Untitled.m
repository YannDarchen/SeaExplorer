 C=NaN(length(explorer2_bydive.pressure),3);
          C(:,1)=explorer2_bydive.p_sorted;
          C(:,2)=explorer2_bydive.temp;
          C(:,3)=explorer2_bydive.s;
          D=sortrows(C,1);
            to_igno = [];
        for i=1:length(D)-1
         if  D(i,1) == D(i+1,1) % delete doublon
           to_igno = [to_igno i+1];
         end
        end
        D(to_igno,:)=[];
        ign = isnan(D(:,1));
        D(ign,:)=[];
        D=real(D);
        