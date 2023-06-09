function flip(x)
    z = zeros(length(x))   
       for i=1:length(x)
           z[i]=x[length(x)+1-i]
       end
   
   return z
   end
   function flipc(x)
   z = zeros(size(x))
   for j=1:size(x)[2]
   
       z[:,j]=flip(x[:,j])
   end
   return z
   end
   function flipr(x)
    z = zeros(size(x))
    for j=1:size(x)[1]
    
        z[j,:]=flip(x[j,:])
    end
    return z
    end