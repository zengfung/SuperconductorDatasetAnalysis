train = readtable('train.csv');
unique = readtable('unique_m.csv');
train = table2array(train);
unique = table2array(unique(:,1:86));

% remove unnecessary columns
idx2keep = [];
for col = 1:86
    if ~isequal(unique(:,col),zeros(21263,1))
        idx2keep = [idx2keep,col];
    end
end
unique = unique(:,idx2keep);

% combine matrices
data = [unique,train(:,1:81)];
% form discrete categories
category = train(:,82);
for i = 1:21263
    if category(i)<50
        category(i) = 1;
    else
        category(i) = 2;
    end
end

%split full data into 3 classes based on categories
c1idx = (category==1);
c2idx = (category==2);
c1 = data(c1idx,:);
c2 = data(c2idx,:);

randidx1 = randperm(14884,floor(0.8*14884));
randidx2 = randperm(6379,floor(0.8*6379));

% training data
train1 = c1(randidx1,:);
train2 = c2(randidx2,:);

% test data
test1 = c1;
test2 = c2;
test1(randidx1,:) = [];
[len_test1,~] = size(test1);
test2(randidx2,:) = [];
[len_test2,~] = size(test2);
test_data = [test1;test2];
test_cat = [1*ones(len_test1,1);2*ones(len_test2,1)];

%% Linear SVM
fprintf("Running linear SVM...\n");
[w_lin,b_lin] = svm(train1,train2,1,2,1e-5,1/21263,'lin');

[test_num,~] = size(test_data);
predlin_val = test_data*w_lin+b_lin;
predlin_cat = zeros(test_num,1);

for i = 1:test_num
    if predlin_val(i) > 0
        predlin_cat(i) = 1;
    elseif predlin_val(i) < 0
        predlin_cat(i) = 2;
    else
        error('predlin_val between -1 and 1');
    end   
end

confmat_lin = confusionmat(test_cat,predlin_cat);
errorrate_lin = (confmat_lin(1,2)+confmat_lin(2,1))/sum(confmat_lin,'all');

%% svm function
function [w,b] = svm(traindata1,traindata2,cat1,cat2,lambda,regparam,type)
    train_data = [traindata1;traindata2];
    [train1_num,~] = size(traindata1);
    [train2_num,~] = size(traindata2);
    train_cat = [cat1*ones(train1_num,1);cat2*ones(train2_num,1)];
    [t_length,w_length] = size(train_data);
    b_length = 1;
    
    D = zeros(t_length); % matrix D (80 x 80) with diagonal values of 1 or -1
    for i = 1:t_length
       if train_cat(i) == cat1
           D(i,i) = 1;
       elseif train_cat(i) == cat2
           D(i,i) = -1;
       else
           %error('invalid data');
           D(i,i) = -1;
       end
    end
    
    e = ones(t_length,1); % vector e of ones
    B = train_data; % matrix A (80 x 15): data matrix for constraints
    m = t_length; % scalar m = number of rows/constraints
    
    H = zeros(t_length+b_length+w_length);
    H(1:w_length,1:w_length) = lambda*eye(w_length);
    f = regparam*[zeros(w_length+b_length,1);e];
    A = -[D*B,D*e,eye(t_length)];
    b = -e;
    LB = [-Inf(w_length+b_length,1);zeros(t_length,1)];
    UB = Inf(w_length+b_length+t_length,1);
    
    if strcmp('lin',type)
        [X,~] = linprog(f,A,b,[],[],LB,UB);
    elseif strcmp('quad',type)
        [X,~] = quadprog(H,f,A,b,[],[],LB,UB);
    else
        error("Invalid program type, select 'lin' or 'quad'");
    end

    w = X(1:w_length);
    b = X(w_length+1);
end