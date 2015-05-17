getMinb <- function(t, out1, regular, npoly=1) {

    if (regular == 'Sparse') {
        dstar <- minb(out1, 2 + npoly); # rough 1D initial value 
# get count matrix 

# 
    
    } else if (regular == 'RegularWithMV') {
        dstar <- minb(out1, 1 + npoly) * 2;
    } else if (regular == 'Dense') {
        dstar = minb(out1, 2 + npoly) * 1.5;
    } 
}


# function [dstar] = getMinb(t,out1,regular, npoly)
# if regular == 0
  # [res] = designPlotCount(t,out1,1,1); %obtain the count matrix based on the design plot 
  # %  fprintf(1,['Time after the count matrix: ' num2str(cputime-t1) '\n']);
  # %  t1 = cputime;

  # dstar = minb(out1,2+npoly);              %rough initial estimate of dstar based on 1-D sense
                                           # %for at least 3 points
  # res(:,[1 end]) = 1;
  # res = res';
  # ids = (res > 0);
  # clear res;
  # b = repmat(out1,length(out1),1)';
 
## WHAT IS THIS?
  # dstar = max(dstar,max(diff(b(ids)))/2);

  # clear b;
# elseif regular == 1
  # dstar = minb(out1,1+npoly)*2;
# else
  # %  t1 = cputime;
  # dstar = minb(out1,2+npoly)*1.5;
# end
# %fprintf(1,['Time after the search of dstar: ' num2str(cputime-t1) '\n']);


# end
