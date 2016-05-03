% Computing variance based upon ABM means and the MLE based upon ABM
% variance

BFR_Young=dlmread('BFR_Raw_Young.txt');

BFR_Old=dlmread('BFR_Raw_Old7.txt');



BFR_Animal=[BFR_Young;BFR_Old];


MLE=0.;

for iprot=1:nprot
    var0=0.;
    for j=1:nnz(BFR_Animal(iprot,:))
        var0=var0+(BFR_Animal(iprot,j)-BFRPeri(iprot,1))^2;
    end
    var0=var0/nnz(BFR_Animal(iprot,:));
    MLE=MLE+nnz(BFR_Animal(iprot,:))/2*(log(var0)+1);
end
MLE=MLE+log(sqrt(2*pi))*nnz(BFR_Animal)% remember, here the negative of the log likelihood is being 
% computed because we are trying to minimize the -ve of the log likelihood
% to approach theoretical MLE using simulated annealing!!!
   
