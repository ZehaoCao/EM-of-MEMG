%Generate network parameters: scripted by Zehao Cao
function [B_out,R,X,n_bus]=LDPF_parameter()
%% Initial power network parameters
mpc = loadcase('case33mg');
%% B matrix generation
R=diag(mpc.branch(:,3));% resistance matrix R
X=diag(mpc.branch(:,4));% reactance matrix X
node = mpc.bus;
n_bus = size(node,1);%number of bus
branch = mpc.branch;
Bra=branch(:,1:2);n=size(Bra,1);
node=zeros(n+1,2);node(:,1)=[Bra(1,1):Bra(1,1)+n];%Record node traversal
bra=zeros(n,2);bra(:,1)=[1:size(bra,1)];%Record node traversal
Road=ones(n+1,n)*NaN;%Record the minimum set of paths
node1=Bra(1,1);
road=zeros(n+1,1);
t=1;
s=1;
road(s)=node1;
while sum(bra(:,2))<n%branches that have not been traversed
    a=find(Bra(:,1)==road(s)|Bra(:,2)==road(s));
    if ~isnan(a)%exist such branches
        b=bra(a,2);
        if sum(b)==size(a,1)%all neighboring branches are traversed
            %backtrack and proceed to next time, record + jump out
            s=s-1;
            c=find(Bra(:,1)==road(s)|Bra(:,2)==road(s));
            if ~isnan(c)
                if t==1
                    Road(:,t)=road;
                    t=t+1;
                else
                mem=ismember(road,Road(:,t-1));
                if sum(mem)~=size(road,1)
                    Road(:,t)=road;
                    t=t+1;
                end
                end
                road=[road(1:s);zeros(size(road,1)-s,1)];
            end
        else
            c=find(bra(a(:,1),2)==0);%finding untraversed branch paths
            c1=Bra(a(c(1)),1:2);
            bra(a(c(1,1)),2)=1;%record that the branch will be traversed soon
            c2=c1(1,:)~=road(s);c2=c1(c2);
            s=s+1;
            road(s)=c2;
            node(c2+1,2)=1;%Record that the node was traversed
        end
    end
end
Road(:,t)=road;
a=find(isnan(Road(1,:)));
Road=Road(:,1:a(1)-1);
B=zeros(n+1,n+1);
for i=2:n+1
    [r1,r2]=find(Road==i);
    B(i,Road(1:r1,r2))=1;
end
B_out = B(2:n_bus,2:n_bus);
end
