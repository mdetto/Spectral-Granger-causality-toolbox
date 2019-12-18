%adjust simultaneously x-lim axes of a multi plot figure m x m
% use this function after multi_spec_GC

function fix_axes(m,xlim,ylim) 

b=0;

      for i=1:m
        for j=1:m
            b=b+1;
          if i~=j  

          subplot(m,m,b)
          set(gca,'ylim',ylim,'xlim',xlim);
          else
          subplot(m,m,b)
          set(gca,'xlim',xlim);              
          end
        end
      end
        