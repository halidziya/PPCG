/*void test1{


cout << "\n";
Vector bet(d);
bet.zero();
bet[0] = 2;
bet[1] = -1;
Normal at(bet, eye(d));
Matrix scat(d, d);
Matrix scot(d, d);
scat.zero();
Vector men(d);
men.zero();
Matrix pis = eye(d);
pis(0, 0) = 2;
pis(0, 1) = 0.2;
pis(1, 0) = 0.2;
IWishart wis(pis, d + 2);
scot.zero();
for (auto i = 0; i < 10000; i++)
{
Vector bit = at.rnd();
men = men + bit;
Vector but = bit - at.mu;
scat = scat + but.outer(but);
scot = scot + wis.rnd();
}
(men / 10000).print();
scat = scat / 10000;
(scot  / 10000).print();
scat.print();
system("pause");



Matrix corr(d,d);
Vector sum(d);
for (int i = 0; i < d; i++)
for (int j = 0; j < d; j++)
corr(j, i) = 1;
Normal norm(v({ 1,1,1,1,1,1,1,1}), eye(d)+corr);
sum.zero();
for (int i = 0; i < 1000; i++)
{
sum = sum + norm.rnd();
}
(sum / 1000).print();



}*/