function Wern=Wern(alpha)

Ze=[1;0];
On=[0;1];

sing=1/sqrt(2)*(kron(Ze,On)-kron(On,Ze));

SING=sing*sing';


Wern=alpha*SING+(1-alpha)*eye(4)/4;