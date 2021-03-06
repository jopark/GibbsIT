function mol=ec_ccm_H(free_net,free_xch,in)

kernel_net=[ -0 -0 -0 -0 -1 -2 2 1 -0 -0 -0 -0 -0 -0 -0 1 -0 -0 -0;
 -0 -0 0.5 0.5 1 1.5 -1.5 -0 1 -0 -0 -0 -0 0.5 -0 -0.5 0.5 -0.5 -0;
 -0 -0 -0 -1 -2 -4 4 -0 -0 1 -0 -0 -0 -0 -0 2 -0 -0 -0;
 -0 -0 -1 -0 -0 -0 -0 -0 -0 -0 1 -0 -0 -0 -0 -0 -1 -1 -0;
 -0 -0 -1 -1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1;
 -0 -0 -0 -0 -1 1 -1 -0 -0 -0 -0 -0 -0 -0 1 -0 -0 -0 -0;
 -1 -1 0.5 0.5 1 -0.5 0.5 -0 -0 -0 -0 1 1 0.5 -0 0.5 0.5 0.5 -0;
 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 -0 -0 0.5 -0.5 -0 -2.5 2.5 -0 -0 -0 -0 -0 1 0.5 -0 1.5 0.5 0.5 -0;
 -0 -0 0.5 -0.5 -0 -0.5 1.5 -0 -0 -0 -0 -0 -0 0.5 -0 0.5 0.5 0.5 -0;
 -0 -0 0.5 -0.5 -0 -0.5 1.5 -0 -0 -0 -0 -0 -0 0.5 -0 0.5 0.5 0.5 -0;
 -0 -0 0.5 -0.5 -0 -0.5 0.5 -0 -0 -0 -0 -0 -0 -0.5 -0 0.5 0.5 0.5 -0;
 -0 -0 1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1 1 -0;
 -0 -0 1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1 1 -0;
 -0 -0 1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1 -0;
 -0 -0 1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1 -0;
 -0 -0 1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1 -0;
 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 -0 -0 -0 1 1 2 -2 -0 -0 -0 -0 -0 -0 -0 -0 -1 -0 -0 -0;
 -0 -0 -0 1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0;
 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 -0 -0 -0 -0 1 2 -2 -0 -0 -0 -0 -0 -0 -0 -0 -1 -0 -0 -0;
 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 -0 -0 -0 -0 -0 2 -2 -0 -0 -0 -0 -0 -0 -0 -0 -1 -0 -0 -0;
 -0 -0 -0 -0 -0 1 -1 0 0 0 0 0 0 0 0 -1 0 0 0;
 -0 -0 -0 -0 -0 1 -1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0;
 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
 -0 -0 -0 -0 -0 -0 1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0;
 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

kernel_xch=[ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];

net=sparse(kernel_net)*free_net;
xch=sparse(kernel_xch)*free_xch;

f=xch+max(0,net);
b=xch-min(0,net);

%%
A1=zeros(100);
A1(1,1)=-b(2)-b(13)-b(18)-b(21)-b(25)-b(29)-f(12)-f(13)-f(19)-f(21)-f(22)-f(23)-f(25)-f(32)-f(10)*f(43)-b(24)*f(44)-b(10)*f(43)-f(24)*f(44);
A1(1,7)=f(13);
A1(1,8)=b(19);
A1(1,14)=b(10)*f(43);
A1(1,24)=b(21);
A1(1,25)=f(21);
A1(1,28)=f(18);
A1(1,29)=f(25);
A1(1,30)=f(29)+b(12);
A1(1,31)=b(13);
A1(1,32)=b(22);
A1(1,33)=b(23)+f(24)*f(44);
A1(1,34)=b(25);
A1(1,43)=f(10)*f(43);
A1(1,67)=b(24)*f(44);
A1(2,2)=-b(3)-b(20)-b(23)-f(33);
A1(2,18)=f(20);
A1(2,21)=f(23);
A1(3,3)=-b(4)-b(14)-f(34);
A1(3,35)=f(14);
A1(4,4)=-b(25)-f(26)-f(27);
A1(4,12)=b(26);
A1(4,13)=b(27);
A1(4,23)=f(25);
A1(5,5)=-b(10)-b(26)-b(28)-f(11)-f(36);
A1(5,9)=f(28);
A1(5,15)=f(26);
A1(5,27)=f(10);
A1(5,36)=b(11);
A1(6,6)=-b(25)-f(26)-f(27);
A1(6,10)=b(27);
A1(6,14)=b(26);
A1(6,33)=f(25);
A1(7,1)=b(13);
A1(7,7)=-b(12)-f(13)-f(29)-f(37);
A1(7,36)=f(12);
A1(7,37)=b(29);
A1(8,1)=f(19)/3+f(22)/3;
A1(8,8)=-b(5)-b(19)-b(22)-f(42);
A1(8,17)=f(19)/3;
A1(8,22)=f(22)/3;
A1(8,24)=f(22)/3;
A1(8,26)=f(19)/3;
A1(9,5)=b(28);
A1(9,9)=-b(27)-b(30)-f(28);
A1(9,16)=f(27);
A1(9,37)=f(30);
A1(10,6)=f(27);
A1(10,10)=-b(27)-b(30)-f(28);
A1(10,14)=b(28);
A1(10,39)=f(30);
A1(11,11)=-b(12)-f(13)-f(29)-f(37);
A1(11,20)=b(13);
A1(11,39)=b(29);
A1(11,40)=f(12);
A1(12,4)=f(26);
A1(12,12)=-b(10)-b(26)-b(28)-f(11)-f(36);
A1(12,13)=f(28);
A1(12,18)=f(10);
A1(12,41)=b(11);
A1(13,4)=f(27);
A1(13,12)=b(28);
A1(13,13)=-b(27)-b(30)-f(28);
A1(13,42)=f(30);
A1(14,1)=f(10)*f(43);
A1(14,6)=f(26);
A1(14,10)=f(28);
A1(14,14)=-b(10)-b(26)-b(28)-f(11)-f(36);
A1(14,40)=b(11);
A1(14,43)=f(10)*(1-f(43));
A1(15,5)=b(26);
A1(15,15)=-b(28)-f(26)-f(29)-f(39);
A1(15,44)=f(28);
A1(15,45)=b(29);
A1(16,9)=b(27);
A1(16,16)=-b(6)-b(24)-f(27)-f(38);
A1(16,23)=f(24);
A1(17,17)=-b(18)-f(19)-f(41);
A1(17,38)=b(19);
A1(17,46)=f(18);
A1(18,2)=b(20);
A1(18,12)=b(10);
A1(18,18)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A1(19,19)=-b(12)-b(13)-b(22)-b(26)-b(27)-f(14)-f(28);
A1(19,47)=f(12);
A1(19,48)=f(13);
A1(19,49)=f(22);
A1(19,50)=f(26)+f(27);
A1(19,51)=b(14);
A1(19,52)=b(28);
A1(20,11)=f(13);
A1(20,20)=-b(12)-b(13)-b(22)-b(26)-b(27)-f(14)-f(28);
A1(20,53)=f(12);
A1(20,54)=f(22);
A1(20,55)=f(26)+f(27);
A1(20,56)=b(14);
A1(20,57)=b(28);
A1(21,2)=b(23);
A1(21,21)=-b(20)-f(21)-f(23);
A1(21,22)=b(21);
A1(21,27)=f(20);
A1(22,8)=b(22);
A1(22,21)=f(21);
A1(22,22)=-b(21)-f(22);
A1(23,4)=b(25);
A1(23,16)=b(24);
A1(23,23)=-b(23)-f(24)-f(25);
A1(23,25)=f(23);
A1(24,1)=f(21);
A1(24,24)=-b(21)-f(22);
A1(24,38)=b(22);
A1(25,1)=b(21);
A1(25,23)=b(23);
A1(25,25)=-b(20)-f(21)-f(23);
A1(25,43)=f(20);
A1(26,26)=-b(18)-f(19)-f(41);
A1(26,32)=b(19);
A1(26,58)=f(18);
A1(27,5)=b(10);
A1(27,21)=b(20);
A1(27,27)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A1(28,1)=b(18);
A1(28,28)=-b(17)-f(18);
A1(28,59)=f(17);
A1(29,1)=b(25);
A1(29,29)=-b(23)-f(24)-f(25);
A1(29,60)=f(23);
A1(29,61)=b(24);
A1(30,1)=f(12)+b(29);
A1(30,30)=-b(12)-f(13)-f(29)-f(37);
A1(30,35)=b(13);
A1(31,1)=f(13);
A1(31,31)=-b(12)-b(13)-b(22)-b(26)-b(27)-f(14)-f(28);
A1(31,62)=f(12);
A1(31,63)=f(22);
A1(31,64)=f(26)+f(27);
A1(31,65)=b(14);
A1(31,66)=b(28);
A1(32,1)=f(19)/3+f(22)/3;
A1(32,17)=f(19)/3;
A1(32,22)=f(22)/3;
A1(32,24)=f(22)/3;
A1(32,26)=f(19)/3;
A1(32,32)=-b(5)-b(19)-b(22)-f(42);
A1(33,1)=f(23)+b(24)*f(44);
A1(33,6)=b(25);
A1(33,33)=-b(23)-f(24)-f(25);
A1(33,67)=b(24)*(1-f(44));
A1(34,1)=f(25);
A1(34,34)=-b(25)-f(26)-f(27);
A1(34,35)=b(26)+b(27);
A1(35,3)=b(14);
A1(35,30)=f(13);
A1(35,34)=f(26)+f(27);
A1(35,35)=-b(12)-b(13)-b(22)-b(26)-b(27)-f(14)-f(28);
A1(35,68)=f(12);
A1(35,69)=f(22);
A1(35,70)=b(28);
A1(36,5)=f(11);
A1(36,7)=b(12);
A1(36,36)=-b(11)-f(12);
A1(37,7)=f(29);
A1(37,9)=b(30);
A1(37,37)=-b(29)-f(30);
A1(38,1)=f(19)/3+f(22)/3;
A1(38,17)=f(19)/3;
A1(38,22)=f(22)/3;
A1(38,24)=f(22)/3;
A1(38,26)=f(19)/3;
A1(38,38)=-b(5)-b(19)-b(22)-f(42);
A1(39,10)=b(30);
A1(39,11)=f(29);
A1(39,39)=-b(29)-f(30);
A1(40,11)=b(12);
A1(40,14)=f(11);
A1(40,40)=-b(11)-f(12);
A1(41,12)=f(11);
A1(41,41)=-b(11)-f(12);
A1(41,48)=b(12);
A1(42,13)=b(30);
A1(42,42)=-b(29)-f(30);
A1(42,48)=f(29);
A1(43,1)=b(10)*f(43);
A1(43,14)=b(10)*(1-f(43));
A1(43,25)=b(20);
A1(43,43)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A1(44,15)=b(28);
A1(44,44)=-b(27)-b(30)-f(28);
A1(44,45)=f(30);
A1(44,67)=f(27);
A1(45,15)=f(29);
A1(45,44)=b(30);
A1(45,45)=-b(29)-f(30);
A1(46,17)=b(18);
A1(46,46)=-b(17)-f(18);
A1(46,71)=f(17);
A1(47,19)=b(12);
A1(47,47)=-b(11)-f(12);
A1(47,52)=f(11);
A1(48,19)=b(13);
A1(48,41)=f(12);
A1(48,42)=b(29);
A1(48,48)=-b(12)-f(13)-f(29)-f(37);
A1(49,19)=b(22);
A1(49,49)=-b(21)-f(22);
A1(49,72)=f(21);
A1(50,19)=b(26)+b(27);
A1(50,50)=-b(25)-f(26)-f(27);
A1(50,73)=f(25);
A1(51,19)=f(14);
A1(51,51)=-b(14)-f(15);
A1(51,74)=b(15);
A1(52,19)=f(28);
A1(52,47)=b(11);
A1(52,52)=-b(10)-b(26)-b(28)-f(11)-f(36);
A1(52,75)=f(10);
A1(52,76)=f(26);
A1(53,20)=b(12);
A1(53,53)=-b(11)-f(12);
A1(53,57)=f(11);
A1(54,20)=b(22);
A1(54,54)=-b(21)-f(22);
A1(54,77)=f(21);
A1(55,20)=b(26)+b(27);
A1(55,55)=-b(25)-f(26)-f(27);
A1(55,78)=f(25);
A1(56,20)=f(14);
A1(56,56)=-b(14)-f(15);
A1(56,79)=b(15);
A1(57,20)=f(28);
A1(57,53)=b(11);
A1(57,57)=-b(10)-b(26)-b(28)-f(11)-f(36);
A1(57,80)=f(10);
A1(57,81)=f(26);
A1(58,26)=b(18);
A1(58,58)=-b(17)-f(18);
A1(58,82)=f(17);
A1(59,28)=b(17);
A1(59,59)=-b(16)-f(17);
A1(59,83)=f(16);
A1(60,29)=b(23);
A1(60,60)=-b(20)-f(21)-f(23);
A1(60,69)=b(21);
A1(60,84)=f(20);
A1(61,29)=f(24);
A1(61,61)=-b(6)-b(24)-f(27)-f(38);
A1(61,85)=b(27);
A1(62,31)=b(12);
A1(62,62)=-b(11)-f(12);
A1(62,66)=f(11);
A1(63,31)=b(22);
A1(63,63)=-b(21)-f(22);
A1(63,86)=f(21);
A1(64,31)=b(26)+b(27);
A1(64,64)=-b(25)-f(26)-f(27);
A1(64,87)=f(25);
A1(65,31)=f(14);
A1(65,65)=-b(14)-f(15);
A1(65,83)=b(15);
A1(66,31)=f(28);
A1(66,62)=b(11);
A1(66,66)=-b(10)-b(26)-b(28)-f(11)-f(36);
A1(66,88)=f(10);
A1(66,89)=f(26);
A1(67,1)=f(24)*f(44);
A1(67,33)=f(24)*(1-f(44));
A1(67,44)=b(27);
A1(67,67)=-b(6)-b(24)-f(27)-f(38);
A1(68,35)=b(12);
A1(68,68)=-b(11)-f(12);
A1(68,70)=f(11);
A1(69,35)=b(22);
A1(69,60)=f(21);
A1(69,69)=-b(21)-f(22);
A1(70,35)=f(28);
A1(70,68)=b(11);
A1(70,70)=-b(10)-b(26)-b(28)-f(11)-f(36);
A1(70,84)=f(10);
A1(70,90)=f(26);
A1(71,46)=b(17);
A1(71,71)=-b(16)-f(17);
A1(71,79)=f(16);
A1(72,49)=b(21);
A1(72,72)=-b(20)-f(21)-f(23);
A1(72,73)=b(23);
A1(72,75)=f(20);
A1(73,50)=b(25);
A1(73,72)=f(23);
A1(73,73)=-b(23)-f(24)-f(25);
A1(73,91)=b(24);
A1(74,51)=f(15);
A1(74,74)=-b(15)-f(16)-f(40);
A1(74,82)=b(16);
A1(75,52)=b(10);
A1(75,72)=b(20);
A1(75,75)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A1(76,52)=b(26);
A1(76,76)=-b(28)-f(26)-f(29)-f(39);
A1(76,92)=f(28);
A1(76,93)=b(29);
A1(77,54)=b(21);
A1(77,77)=-b(20)-f(21)-f(23);
A1(77,78)=b(23);
A1(77,80)=f(20);
A1(78,55)=b(25);
A1(78,77)=f(23);
A1(78,78)=-b(23)-f(24)-f(25);
A1(78,94)=b(24);
A1(79,56)=f(15);
A1(79,71)=b(16);
A1(79,79)=-b(15)-f(16)-f(40);
A1(80,57)=b(10);
A1(80,77)=b(20);
A1(80,80)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A1(81,57)=b(26);
A1(81,81)=-b(28)-f(26)-f(29)-f(39);
A1(81,95)=f(28);
A1(81,96)=b(29);
A1(82,58)=b(17);
A1(82,74)=f(16);
A1(82,82)=-b(16)-f(17);
A1(83,59)=b(16);
A1(83,65)=f(15);
A1(83,83)=-b(15)-f(16)-f(40);
A1(84,60)=b(20);
A1(84,70)=b(10);
A1(84,84)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A1(85,61)=f(27);
A1(85,85)=-b(27)-b(30)-f(28);
A1(85,90)=b(28);
A1(85,97)=f(30);
A1(86,63)=b(21);
A1(86,86)=-b(20)-f(21)-f(23);
A1(86,87)=b(23);
A1(86,88)=f(20);
A1(87,64)=b(25);
A1(87,86)=f(23);
A1(87,87)=-b(23)-f(24)-f(25);
A1(87,98)=b(24);
A1(88,66)=b(10);
A1(88,86)=b(20);
A1(88,88)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A1(89,66)=b(26);
A1(89,89)=-b(28)-f(26)-f(29)-f(39);
A1(89,99)=f(28);
A1(89,100)=b(29);
A1(90,70)=b(26);
A1(90,85)=f(28);
A1(90,90)=-b(28)-f(26)-f(29)-f(39);
A1(90,97)=b(29);
A1(91,73)=f(24);
A1(91,91)=-b(6)-b(24)-f(27)-f(38);
A1(91,92)=b(27);
A1(92,76)=b(28);
A1(92,91)=f(27);
A1(92,92)=-b(27)-b(30)-f(28);
A1(92,93)=f(30);
A1(93,76)=f(29);
A1(93,92)=b(30);
A1(93,93)=-b(29)-f(30);
A1(94,78)=f(24);
A1(94,94)=-b(6)-b(24)-f(27)-f(38);
A1(94,95)=b(27);
A1(95,81)=b(28);
A1(95,94)=f(27);
A1(95,95)=-b(27)-b(30)-f(28);
A1(95,96)=f(30);
A1(96,81)=f(29);
A1(96,95)=b(30);
A1(96,96)=-b(29)-f(30);
A1(97,85)=b(30);
A1(97,90)=f(29);
A1(97,97)=-b(29)-f(30);
A1(98,87)=f(24);
A1(98,98)=-b(6)-b(24)-f(27)-f(38);
A1(98,99)=b(27);
A1(99,89)=b(28);
A1(99,98)=f(27);
A1(99,99)=-b(27)-b(30)-f(28);
A1(99,100)=f(30);
A1(100,89)=f(29);
A1(100,99)=b(30);
A1(100,100)=-b(29)-f(30);

B1=zeros(100,33);
B1(1,1)=-f(2);
B1(2,2)=-f(3);
B1(3,3)=-f(4);
B1(8,4)=-f(5);
B1(16,5)=-f(6);
B1(18,6)=-f(7);
B1(18,7)=-f(8);
B1(18,8)=-f(9);
B1(27,9)=-f(7);
B1(27,10)=-f(8);
B1(27,11)=-f(9);
B1(32,12)=-f(5);
B1(38,13)=-f(5);
B1(43,14)=-f(7);
B1(43,15)=-f(8);
B1(43,16)=-f(9);
B1(61,17)=-f(6);
B1(67,18)=-f(6);
B1(75,19)=-f(7);
B1(75,20)=-f(8);
B1(75,21)=-f(9);
B1(80,22)=-f(7);
B1(80,23)=-f(8);
B1(80,24)=-f(9);
B1(84,25)=-f(7);
B1(84,26)=-f(8);
B1(84,27)=-f(9);
B1(88,28)=-f(7);
B1(88,29)=-f(8);
B1(88,30)=-f(9);
B1(91,31)=-f(6);
B1(94,32)=-f(6);
B1(98,33)=-f(6);

y1=[ in.H_IN_1; in.NADPH_IN_1; in.NADH_IN_1; in.PYR_IN_6; in.RIBOSE_IN_6; in.GLC_IN_7; in.GLYCOGEN_IN1_7; in.GLYCOGEN_IN2_7; in.GLC_IN_9; in.GLYCOGEN_IN1_9; in.GLYCOGEN_IN2_9; in.PYR_IN_4; in.PYR_IN_5; in.GLC_IN_8; in.GLYCOGEN_IN1_8; in.GLYCOGEN_IN2_8; in.RIBOSE_IN_8; in.RIBOSE_IN_7; in.GLC_IN_12; in.GLYCOGEN_IN1_12; in.GLYCOGEN_IN2_12; in.GLC_IN_13; in.GLYCOGEN_IN1_13; in.GLYCOGEN_IN2_13; in.GLC_IN_10; in.GLYCOGEN_IN1_10; in.GLYCOGEN_IN2_10; in.GLC_IN_11; in.GLYCOGEN_IN1_11; in.GLYCOGEN_IN2_11; in.RIBOSE_IN_10; in.RIBOSE_IN_11; in.RIBOSE_IN_9];
x1=sparse(A1)\(B1*sparse(y1));

%%
A2=zeros(32);
A2(1,1)=-b(18)-f(19)-f(41);
A2(1,14)=f(18);
A2(1,15)=b(19);
A2(2,2)=-b(25)-f(26)-f(27);
A2(2,4)=b(26);
A2(2,5)=b(27);
A2(2,16)=f(25);
A2(3,1)=f(19)/3;
A2(3,3)=-b(5)-b(19)-b(22)-f(42);
A2(3,6)=f(22)/3;
A2(4,2)=f(26);
A2(4,4)=-b(10)-b(26)-b(28)-f(11)-f(36);
A2(4,5)=f(28);
A2(4,9)=b(11);
A2(4,18)=f(10)*(1-f(43));
A2(5,2)=f(27);
A2(5,4)=b(28);
A2(5,5)=-b(27)-b(30)-f(28);
A2(5,10)=f(30);
A2(6,3)=b(22);
A2(6,6)=-b(21)-f(22);
A2(7,7)=-b(12)-b(13)-b(22)-b(26)-b(27)-f(14)-f(28);
A2(7,8)=f(13);
A2(7,11)=b(14);
A2(7,13)=f(26)+f(27);
A2(7,19)=f(12);
A2(7,20)=f(22);
A2(7,21)=b(28);
A2(8,7)=b(13);
A2(8,8)=-b(12)-f(13)-f(29)-f(37);
A2(8,9)=f(12);
A2(8,10)=b(29);
A2(9,4)=f(11);
A2(9,8)=b(12);
A2(9,9)=-b(11)-f(12);
A2(10,5)=b(30);
A2(10,8)=f(29);
A2(10,10)=-b(29)-f(30);
A2(11,7)=f(14);
A2(11,11)=-b(14)-f(15);
A2(11,22)=b(15);
A2(12,12)=-b(23)-f(24)-f(25);
A2(12,13)=b(25);
A2(12,23)=f(23);
A2(12,24)=b(24);
A2(13,7)=b(26)+b(27);
A2(13,12)=f(25);
A2(13,13)=-b(25)-f(26)-f(27);
A2(14,1)=b(18);
A2(14,14)=-b(17)-f(18);
A2(14,25)=f(17);
A2(15,1)=f(19)/3;
A2(15,6)=f(22)/3;
A2(15,15)=-b(5)-b(19)-b(22)-f(42);
A2(16,2)=b(25);
A2(16,16)=-b(23)-f(24)-f(25);
A2(16,26)=b(24)*(1-f(44));
A2(17,1)=f(19)/3;
A2(17,6)=f(22)/3;
A2(17,17)=-b(5)-b(19)-b(22)-f(42);
A2(18,4)=b(10)*(1-f(43));
A2(18,18)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A2(19,7)=b(12);
A2(19,19)=-b(11)-f(12);
A2(19,21)=f(11);
A2(20,7)=b(22);
A2(20,20)=-b(21)-f(22);
A2(20,23)=f(21);
A2(21,7)=f(28);
A2(21,19)=b(11);
A2(21,21)=-b(10)-b(26)-b(28)-f(11)-f(36);
A2(21,27)=f(10);
A2(21,28)=f(26);
A2(22,11)=f(15);
A2(22,22)=-b(15)-f(16)-f(40);
A2(22,25)=b(16);
A2(23,12)=b(23);
A2(23,20)=b(21);
A2(23,23)=-b(20)-f(21)-f(23);
A2(23,27)=f(20);
A2(24,12)=f(24);
A2(24,24)=-b(6)-b(24)-f(27)-f(38);
A2(24,29)=b(27);
A2(25,14)=b(17);
A2(25,22)=f(16);
A2(25,25)=-b(16)-f(17);
A2(26,16)=f(24)*(1-f(44));
A2(26,26)=-b(6)-b(24)-f(27)-f(38);
A2(26,30)=b(27);
A2(27,21)=b(10);
A2(27,23)=b(20);
A2(27,27)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A2(28,21)=b(26);
A2(28,28)=-b(28)-f(26)-f(29)-f(39);
A2(28,29)=f(28);
A2(28,31)=b(29);
A2(29,24)=f(27);
A2(29,28)=b(28);
A2(29,29)=-b(27)-b(30)-f(28);
A2(29,31)=f(30);
A2(30,26)=f(27);
A2(30,30)=-b(27)-b(30)-f(28);
A2(30,32)=f(30);
A2(31,28)=f(29);
A2(31,29)=b(30);
A2(31,31)=-b(29)-f(30);
A2(32,30)=b(30);
A2(32,32)=-b(29)-f(30);

B2=zeros(32,32);
B2(3,1)=-f(19)/3;
B2(3,2)=-f(19)/3;
B2(3,3)=-f(22)/3;
B2(3,4)=-f(22)/3;
B2(3,17)=-f(5);
B2(4,5)=-f(10)*f(43);
B2(6,6)=-f(21);
B2(15,7)=-f(19)/3;
B2(15,8)=-f(19)/3;
B2(15,9)=-f(22)/3;
B2(15,10)=-f(22)/3;
B2(15,20)=-f(5);
B2(16,11)=-f(23);
B2(16,12)=-b(24)*f(44);
B2(17,13)=-f(19)/3;
B2(17,14)=-f(19)/3;
B2(17,15)=-f(22)/3;
B2(17,16)=-f(22)/3;
B2(17,22)=-f(5);
B2(18,18)=-b(10)*f(43);
B2(18,19)=-b(20);
B2(18,23)=-f(7);
B2(18,24)=-f(8);
B2(18,25)=-f(9);
B2(24,26)=-f(6);
B2(26,21)=-f(24)*f(44);
B2(26,27)=-f(6);
B2(27,29)=-f(7);
B2(27,30)=-f(8);
B2(27,31)=-f(9);
B2(30,28)=-b(28);
B2(32,32)=-f(29);

y2=[ cauchy(x1(1,:),x1(17,:)); cauchy(x1(1,:),x1(26,:)); cauchy(x1(1,:),x1(24,:)); cauchy(x1(1,:),x1(22,:)); cauchy(x1(18,:),x1(1,:)); cauchy(x1(1,:),x1(21,:)); cauchy(x1(1,:),x1(17,:)); cauchy(x1(1,:),x1(26,:)); cauchy(x1(1,:),x1(24,:)); cauchy(x1(1,:),x1(22,:)); cauchy(x1(1,:),x1(25,:)); cauchy(x1(16,:),x1(1,:)); cauchy(x1(1,:),x1(17,:)); cauchy(x1(1,:),x1(26,:)); cauchy(x1(1,:),x1(24,:)); cauchy(x1(1,:),x1(22,:)); in.PYR_IN_5_6; cauchy(x1(12,:),x1(1,:)); cauchy(x1(2,:),x1(25,:)); in.PYR_IN_4_5; cauchy(x1(23,:),x1(1,:)); in.PYR_IN_4_6; in.GLC_IN_7_8; in.GLYCOGEN_IN1_7_8; in.GLYCOGEN_IN2_7_8; in.RIBOSE_IN_10_11; in.RIBOSE_IN_6_7; cauchy(x1(5,:),x1(15,:)); in.GLC_IN_12_13; in.GLYCOGEN_IN1_12_13; in.GLYCOGEN_IN2_12_13; cauchy(x1(7,:),x1(15,:))];
x2=sparse(A2)\(B2*sparse(y2));

%%
A3=zeros(6);
A3(1,1)=-b(27)-b(30)-f(28);
A3(1,4)=b(28);
A3(1,5)=f(30);
A3(2,2)=-b(12)-f(13)-f(29)-f(37);
A3(2,5)=b(29);
A3(2,6)=f(12);
A3(3,3)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A3(3,4)=b(10);
A3(4,1)=f(28);
A3(4,3)=f(10);
A3(4,4)=-b(10)-b(26)-b(28)-f(11)-f(36);
A3(4,6)=b(11);
A3(5,1)=b(30);
A3(5,2)=f(29);
A3(5,5)=-b(29)-f(30);
A3(6,2)=b(12);
A3(6,4)=f(11);
A3(6,6)=-b(11)-f(12);

B3=zeros(6,7);
B3(1,1)=-f(27);
B3(2,2)=-b(13);
B3(3,3)=-b(20);
B3(3,5)=-f(7);
B3(3,6)=-f(8);
B3(3,7)=-f(9);
B3(4,4)=-f(26);

y3=[ cauchy(x1(4,:),x1(16,:)); cauchy(x1(1,:),x1(19,:)); cauchy(x1(2,:),x1(21,:)); cauchy(x1(4,:),x1(15,:)); in.GLC_IN_7_9; in.GLYCOGEN_IN1_7_9; in.GLYCOGEN_IN2_7_9];
x3=A3\(B3*sparse(y3));

%%
A4=zeros(7);
A4(1,1)=-b(27)-b(30)-f(28);
A4(1,4)=f(30);
A4(1,5)=b(28);
A4(2,2)=-b(12)-f(13)-f(29)-f(37);
A4(2,4)=b(29);
A4(2,6)=f(12);
A4(3,3)=-b(20)-f(21)-f(23);
A4(3,7)=f(20);
A4(4,1)=b(30);
A4(4,2)=f(29);
A4(4,4)=-b(29)-f(30);
A4(5,1)=f(28);
A4(5,5)=-b(10)-b(26)-b(28)-f(11)-f(36);
A4(5,6)=b(11);
A4(5,7)=f(10)*(1-f(43));
A4(6,2)=b(12);
A4(6,5)=f(11);
A4(6,6)=-b(11)-f(12);
A4(7,3)=b(20);
A4(7,5)=b(10)*(1-f(43));
A4(7,7)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);

B4=zeros(7,10);
B4(1,1)=-f(27);
B4(2,2)=-b(13);
B4(3,3)=-b(21);
B4(3,4)=-b(23);
B4(5,5)=-f(10)*f(43);
B4(5,6)=-f(26);
B4(7,7)=-b(10)*f(43);
B4(7,8)=-f(7);
B4(7,9)=-f(8);
B4(7,10)=-f(9);

y4=[ cauchy(x1(6,:),x1(16,:)); cauchy(x1(1,:),x1(20,:)); cauchy(x1(1,:),x1(22,:)); cauchy(x1(2,:),x1(23,:)); cauchy(x1(27,:),x1(1,:)); cauchy(x1(6,:),x1(15,:)); cauchy(x1(5,:),x1(1,:)); in.GLC_IN_8_9; in.GLYCOGEN_IN1_8_9; in.GLYCOGEN_IN2_8_9];
x4=A4\(B4*sparse(y4));

%%
A5=zeros(16);
A5(1,1)=-b(14)-f(15);
A5(1,2)=b(15);
A5(1,5)=f(14);
A5(2,1)=f(15);
A5(2,2)=-b(15)-f(16)-f(40);
A5(2,3)=b(16);
A5(3,2)=f(16);
A5(3,3)=-b(16)-f(17);
A5(3,4)=b(17);
A5(4,3)=f(17);
A5(4,4)=-b(17)-f(18);
A5(5,1)=b(14);
A5(5,5)=-b(12)-b(13)-b(22)-b(26)-b(27)-f(14)-f(28);
A5(5,6)=f(26)+f(27);
A5(5,7)=f(22);
A5(5,9)=b(28);
A5(5,14)=f(12);
A5(6,5)=b(26)+b(27);
A5(6,6)=-b(25)-f(26)-f(27);
A5(6,8)=f(25);
A5(7,5)=b(22);
A5(7,7)=-b(21)-f(22);
A5(7,11)=f(21);
A5(8,6)=b(25);
A5(8,8)=-b(23)-f(24)-f(25);
A5(8,11)=f(23);
A5(8,12)=b(24);
A5(9,5)=f(28);
A5(9,9)=-b(10)-b(26)-b(28)-f(11)-f(36);
A5(9,10)=f(26);
A5(9,13)=f(10);
A5(9,14)=b(11);
A5(10,9)=b(26);
A5(10,10)=-b(28)-f(26)-f(29)-f(39);
A5(10,15)=f(28);
A5(10,16)=b(29);
A5(11,7)=b(21);
A5(11,8)=b(23);
A5(11,11)=-b(20)-f(21)-f(23);
A5(11,13)=f(20);
A5(12,8)=f(24);
A5(12,12)=-b(6)-b(24)-f(27)-f(38);
A5(12,15)=b(27);
A5(13,9)=b(10);
A5(13,11)=b(20);
A5(13,13)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A5(14,5)=b(12);
A5(14,9)=f(11);
A5(14,14)=-b(11)-f(12);
A5(15,10)=b(28);
A5(15,12)=f(27);
A5(15,15)=-b(27)-b(30)-f(28);
A5(15,16)=f(30);
A5(16,10)=f(29);
A5(16,15)=b(30);
A5(16,16)=-b(29)-f(30);

B5=zeros(16,6);
B5(4,1)=-b(18);
B5(5,2)=-f(13);
B5(12,3)=-f(6);
B5(13,4)=-f(7);
B5(13,5)=-f(8);
B5(13,6)=-f(9);

y5=[ cauchy(x1(1,:),x2(1,:)); cauchy(x1(1,:),x2(8,:)); in.RIBOSE_IN_9_10_11; in.GLC_IN_11_12_13; in.GLYCOGEN_IN1_11_12_13; in.GLYCOGEN_IN2_11_12_13];
x5=sparse(A5)\(B5*sparse(y5));

%%
A6=zeros(1);
A6(1,1)=-b(5)-b(19)-b(22)-f(42);

B6=zeros(1,3);
B6(1,1)=-f(19);
B6(1,2)=-f(22);
B6(1,3)=-f(5);

y6=[ cauchy(x1(1,:),x2(1,:)); cauchy(x1(1,:),x2(6,:)); in.PYR_IN_4_5_6];
x6=A6\(B6*sparse(y6));

%%
A7=zeros(6);
A7(1,1)=-b(10)-b(26)-b(28)-f(11)-f(36);
A7(1,3)=f(28);
A7(1,4)=b(11);
A7(1,6)=f(10)*(1-f(43));
A7(2,2)=-b(12)-f(13)-f(29)-f(37);
A7(2,4)=f(12);
A7(2,5)=b(29);
A7(3,1)=b(28);
A7(3,3)=-b(27)-b(30)-f(28);
A7(3,5)=f(30);
A7(4,1)=f(11);
A7(4,2)=b(12);
A7(4,4)=-b(11)-f(12);
A7(5,2)=f(29);
A7(5,3)=b(30);
A7(5,5)=-b(29)-f(30);
A7(6,1)=b(10)*(1-f(43));
A7(6,6)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);

B7=zeros(6,9);
B7(1,1)=-f(10)*f(43);
B7(1,2)=-f(26);
B7(2,3)=-b(13);
B7(3,4)=-f(27);
B7(6,5)=-b(10)*f(43);
B7(6,6)=-b(20);
B7(6,7)=-f(7);
B7(6,8)=-f(8);
B7(6,9)=-f(9);

y7=[ cauchy(x3(3,:),x1(1,:)); cauchy(x2(2,:),x1(15,:)); cauchy(x1(1,:),x2(7,:)); cauchy(x2(2,:),x1(16,:)); cauchy(x3(4,:),x1(1,:)); cauchy(x1(2,:),x4(3,:)); in.GLC_IN_7_8_9; in.GLYCOGEN_IN1_7_8_9; in.GLYCOGEN_IN2_7_8_9];
x7=A7\(B7*sparse(y7));

%%
A8=zeros(13);
A8(1,1)=-b(12)-f(13)-f(29)-f(37);
A8(1,2)=b(13);
A8(2,1)=f(13);
A8(2,2)=-b(12)-b(13)-b(22)-b(26)-b(27)-f(14)-f(28);
A8(2,3)=f(12);
A8(2,4)=f(22);
A8(2,5)=f(26)+f(27);
A8(2,6)=b(28);
A8(3,2)=b(12);
A8(3,3)=-b(11)-f(12);
A8(3,6)=f(11);
A8(4,2)=b(22);
A8(4,4)=-b(21)-f(22);
A8(4,7)=f(21);
A8(5,2)=b(26)+b(27);
A8(5,5)=-b(25)-f(26)-f(27);
A8(6,2)=f(28);
A8(6,3)=b(11);
A8(6,6)=-b(10)-b(26)-b(28)-f(11)-f(36);
A8(6,8)=f(10);
A8(6,9)=f(26);
A8(7,4)=b(21);
A8(7,7)=-b(20)-f(21)-f(23);
A8(7,8)=f(20);
A8(7,10)=b(23);
A8(8,6)=b(10);
A8(8,7)=b(20);
A8(8,8)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A8(9,6)=b(26);
A8(9,9)=-b(28)-f(26)-f(29)-f(39);
A8(9,11)=f(28);
A8(9,12)=b(29);
A8(10,7)=f(23);
A8(10,10)=-b(23)-f(24)-f(25);
A8(10,13)=b(24);
A8(11,9)=b(28);
A8(11,11)=-b(27)-b(30)-f(28);
A8(11,12)=f(30);
A8(11,13)=f(27);
A8(12,9)=f(29);
A8(12,11)=b(30);
A8(12,12)=-b(29)-f(30);
A8(13,10)=f(24);
A8(13,11)=b(27);
A8(13,13)=-b(6)-b(24)-f(27)-f(38);

B8=zeros(13,9);
B8(1,1)=-f(12);
B8(1,2)=-b(29);
B8(2,3)=-b(14);
B8(5,4)=-f(25);
B8(8,6)=-f(7);
B8(8,7)=-f(8);
B8(8,8)=-f(9);
B8(10,5)=-b(25);
B8(13,9)=-f(6);

y8=[ cauchy(x1(1,:),x2(9,:)); cauchy(x1(1,:),x2(10,:)); cauchy(x1(3,:),x2(11,:)); cauchy(x1(1,:),x2(12,:)); cauchy(x1(1,:),x2(13,:)); in.GLC_IN_10_12_13; in.GLYCOGEN_IN1_10_12_13; in.GLYCOGEN_IN2_10_12_13; in.RIBOSE_IN_8_10_11];
x8=sparse(A8)\(B8*sparse(y8));

%%
A9=zeros(12);
A9(1,1)=-b(12)-b(13)-b(22)-b(26)-b(27)-f(14)-f(28);
A9(1,3)=f(22);
A9(1,4)=b(28);
A9(1,9)=f(12);
A9(1,10)=f(26)+f(27);
A9(2,2)=-b(23)-f(24)-f(25);
A9(2,7)=f(23);
A9(2,8)=b(24);
A9(3,1)=b(22);
A9(3,3)=-b(21)-f(22);
A9(3,7)=f(21);
A9(4,1)=f(28);
A9(4,4)=-b(10)-b(26)-b(28)-f(11)-f(36);
A9(4,5)=f(26);
A9(4,6)=f(10);
A9(4,9)=b(11);
A9(5,4)=b(26);
A9(5,5)=-b(28)-f(26)-f(29)-f(39);
A9(5,11)=f(28);
A9(5,12)=b(29);
A9(6,4)=b(10);
A9(6,6)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A9(6,7)=b(20);
A9(7,2)=b(23);
A9(7,3)=b(21);
A9(7,6)=f(20);
A9(7,7)=-b(20)-f(21)-f(23);
A9(8,2)=f(24);
A9(8,8)=-b(6)-b(24)-f(27)-f(38);
A9(8,11)=b(27);
A9(9,1)=b(12);
A9(9,4)=f(11);
A9(9,9)=-b(11)-f(12);
A9(10,1)=b(26)+b(27);
A9(10,10)=-b(25)-f(26)-f(27);
A9(11,5)=b(28);
A9(11,8)=f(27);
A9(11,11)=-b(27)-b(30)-f(28);
A9(11,12)=f(30);
A9(12,5)=f(29);
A9(12,11)=b(30);
A9(12,12)=-b(29)-f(30);

B9=zeros(12,8);
B9(1,1)=-f(13);
B9(1,2)=-b(14);
B9(2,3)=-b(25);
B9(6,5)=-f(7);
B9(6,6)=-f(8);
B9(6,7)=-f(9);
B9(8,8)=-f(6);
B9(10,4)=-f(25);

y9=[ cauchy(x1(1,:),x8(1,:)); cauchy(x1(3,:),x5(1,:)); cauchy(x1(1,:),x5(6,:)); cauchy(x1(1,:),x5(8,:)); in.GLC_IN_10_11_12_13; in.GLYCOGEN_IN1_10_11_12_13; in.GLYCOGEN_IN2_10_11_12_13; in.RIBOSE_IN_8_9_10_11];
x9=sparse(A9)\(B9*sparse(y9));

%%
A10=zeros(1);
A10(1,1)=-b(12)-f(13)-f(29)-f(37);

B10=zeros(1,3);
B10(1,1)=-f(12);
B10(1,2)=-b(13);
B10(1,3)=-b(29);

y10=[ cauchy(x1(1,:),x7(4,:)); cauchy(x1(1,:),x8(2,:)); cauchy(x1(1,:),x7(5,:))];
x10=A10\(B10*sparse(y10));

%%
A11=zeros(9);
A11(1,1)=-b(25)-f(26)-f(27);
A11(1,4)=f(25);
A11(2,2)=-b(20)-f(21)-f(23);
A11(2,4)=b(23);
A11(2,5)=f(20);
A11(3,3)=-b(6)-b(24)-f(27)-f(38);
A11(3,4)=f(24);
A11(3,6)=b(27);
A11(4,1)=b(25);
A11(4,2)=f(23);
A11(4,3)=b(24);
A11(4,4)=-b(23)-f(24)-f(25);
A11(5,2)=b(20);
A11(5,5)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A11(5,7)=b(10)*(1-f(43));
A11(6,3)=f(27);
A11(6,6)=-b(27)-b(30)-f(28);
A11(6,8)=f(30);
A11(7,5)=f(10)*(1-f(43));
A11(7,7)=-b(10)-b(26)-b(28)-f(11)-f(36);
A11(7,9)=b(11);
A11(8,6)=b(30);
A11(8,8)=-b(29)-f(30);
A11(9,7)=f(11);
A11(9,9)=-b(11)-f(12);

B11=zeros(9,14);
B11(1,1)=-b(26);
B11(1,2)=-b(27);
B11(2,3)=-b(21);
B11(3,5)=-f(6);
B11(5,4)=-b(10)*f(43);
B11(5,7)=-f(7);
B11(5,8)=-f(8);
B11(5,9)=-f(9);
B11(6,6)=-b(28);
B11(7,10)=-f(10)*f(43);
B11(7,11)=-f(26);
B11(7,12)=-f(28);
B11(8,13)=-f(29);
B11(9,14)=-b(12);

y11=[ cauchy(x1(12,:),x5(5,:)); cauchy(x1(13,:),x5(5,:)); cauchy(x1(1,:),x5(7,:)); cauchy(x5(9,:),x1(1,:)); in.RIBOSE_IN_6_9_10_11; cauchy(x1(5,:),x5(10,:)); in.GLC_IN_8_11_12_13; in.GLYCOGEN_IN1_8_11_12_13; in.GLYCOGEN_IN2_8_11_12_13; cauchy(x5(13,:),x1(1,:)); cauchy(x1(6,:),x5(10,:)); cauchy(x1(10,:),x5(5,:)); cauchy(x1(7,:),x5(10,:)); cauchy(x1(11,:),x5(5,:))];
x11=A11\(B11*sparse(y11));

%%
A12=zeros(11);
A12(1,1)=-b(28)-f(26)-f(29)-f(39);
A12(1,3)=f(28);
A12(1,4)=b(26);
A12(1,5)=b(29);
A12(2,2)=-b(25)-f(26)-f(27);
A12(2,6)=f(25);
A12(3,1)=b(28);
A12(3,3)=-b(27)-b(30)-f(28);
A12(3,5)=f(30);
A12(3,7)=f(27);
A12(4,1)=f(26);
A12(4,4)=-b(10)-b(26)-b(28)-f(11)-f(36);
A12(4,8)=f(10);
A12(4,9)=b(11);
A12(5,1)=f(29);
A12(5,3)=b(30);
A12(5,5)=-b(29)-f(30);
A12(6,2)=b(25);
A12(6,6)=-b(23)-f(24)-f(25);
A12(6,7)=b(24)*(1-f(44));
A12(7,3)=b(27);
A12(7,6)=f(24)*(1-f(44));
A12(7,7)=-b(6)-b(24)-f(27)-f(38);
A12(8,4)=b(10);
A12(8,8)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A12(8,10)=b(20);
A12(9,4)=f(11);
A12(9,9)=-b(11)-f(12);
A12(10,8)=f(20);
A12(10,10)=-b(20)-f(21)-f(23);
A12(10,11)=b(21);
A12(11,10)=f(21);
A12(11,11)=-b(21)-f(22);

B12=zeros(11,13);
B12(2,1)=-b(26);
B12(2,2)=-b(27);
B12(4,3)=-f(28);
B12(6,4)=-f(23);
B12(6,5)=-b(24)*f(44);
B12(7,6)=-f(24)*f(44);
B12(7,8)=-f(6);
B12(8,9)=-f(7);
B12(8,10)=-f(8);
B12(8,11)=-f(9);
B12(9,7)=-b(12);
B12(10,12)=-b(23);
B12(11,13)=-b(22);

y12=[ cauchy(x1(14,:),x5(5,:)); cauchy(x1(10,:),x5(5,:)); cauchy(x1(9,:),x5(5,:)); cauchy(x1(1,:),x5(11,:)); cauchy(x5(12,:),x1(1,:)); cauchy(x5(8,:),x1(1,:)); cauchy(x1(7,:),x5(5,:)); in.RIBOSE_IN_7_9_10_11; in.GLC_IN_9_11_12_13; in.GLYCOGEN_IN1_9_11_12_13; in.GLYCOGEN_IN2_9_11_12_13; cauchy(x1(2,:),x5(8,:)); cauchy(x1(8,:),x5(5,:))];
x12=sparse(A12)\(B12*sparse(y12));

%%
A13=zeros(10);
A13(1,1)=-b(28)-f(26)-f(29)-f(39);
A13(1,4)=b(26);
A13(1,6)=f(28);
A13(1,7)=b(29);
A13(2,2)=-b(21)-f(22);
A13(2,3)=f(21);
A13(3,2)=b(21);
A13(3,3)=-b(20)-f(21)-f(23);
A13(3,5)=f(20);
A13(4,1)=f(26);
A13(4,4)=-b(10)-b(26)-b(28)-f(11)-f(36);
A13(4,5)=f(10);
A13(4,8)=b(11);
A13(5,3)=b(20);
A13(5,4)=b(10);
A13(5,5)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A13(6,1)=b(28);
A13(6,6)=-b(27)-b(30)-f(28);
A13(6,7)=f(30);
A13(6,9)=f(27);
A13(7,1)=f(29);
A13(7,6)=b(30);
A13(7,7)=-b(29)-f(30);
A13(8,4)=f(11);
A13(8,8)=-b(11)-f(12);
A13(9,6)=b(27);
A13(9,9)=-b(6)-b(24)-f(27)-f(38);
A13(9,10)=f(24)*(1-f(44));
A13(10,9)=b(24)*(1-f(44));
A13(10,10)=-b(23)-f(24)-f(25);

B13=zeros(10,12);
B13(2,1)=-b(22);
B13(3,2)=-b(23);
B13(4,3)=-f(28);
B13(5,5)=-f(7);
B13(5,6)=-f(8);
B13(5,7)=-f(9);
B13(8,4)=-b(12);
B13(9,8)=-f(24)*f(44);
B13(9,9)=-f(6);
B13(10,10)=-f(23);
B13(10,11)=-b(24)*f(44);
B13(10,12)=-b(25);

y13=[ cauchy(x1(8,:),x9(1,:)); cauchy(x1(2,:),x9(2,:)); cauchy(x1(9,:),x9(1,:)); cauchy(x1(7,:),x9(1,:)); in.GLC_IN_9_10_11_12_13; in.GLYCOGEN_IN1_9_10_11_12_13; in.GLYCOGEN_IN2_9_10_11_12_13; cauchy(x9(2,:),x1(1,:)); in.RIBOSE_IN_7_8_9_10_11; cauchy(x1(1,:),x9(7,:)); cauchy(x9(8,:),x1(1,:)); cauchy(x1(1,:),x12(2,:))];
x13=A13\(B13*sparse(y13));

%%
A14=zeros(8);
A14(1,1)=-b(23)-f(24)-f(25);
A14(1,2)=f(23);
A14(1,3)=b(24);
A14(2,1)=b(23);
A14(2,2)=-b(20)-f(21)-f(23);
A14(2,4)=f(20);
A14(3,1)=f(24);
A14(3,3)=-b(6)-b(24)-f(27)-f(38);
A14(3,5)=b(27);
A14(4,2)=b(20);
A14(4,4)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A14(4,6)=b(10)*(1-f(43));
A14(5,3)=f(27);
A14(5,5)=-b(27)-b(30)-f(28);
A14(5,7)=f(30);
A14(6,4)=f(10)*(1-f(43));
A14(6,6)=-b(10)-b(26)-b(28)-f(11)-f(36);
A14(6,8)=b(11);
A14(7,5)=b(30);
A14(7,7)=-b(29)-f(30);
A14(8,6)=f(11);
A14(8,8)=-b(11)-f(12);

B14=zeros(8,13);
B14(1,1)=-b(25);
B14(2,2)=-b(21);
B14(3,4)=-f(6);
B14(4,3)=-b(10)*f(43);
B14(4,6)=-f(7);
B14(4,7)=-f(8);
B14(4,8)=-f(9);
B14(5,5)=-b(28);
B14(6,9)=-f(10)*f(43);
B14(6,10)=-f(26);
B14(6,11)=-f(28);
B14(7,12)=-f(29);
B14(8,13)=-b(12);

y14=[ cauchy(x1(1,:),x11(1,:)); cauchy(x1(1,:),x9(3,:)); cauchy(x9(4,:),x1(1,:)); in.RIBOSE_IN_6_8_9_10_11; cauchy(x1(5,:),x9(5,:)); in.GLC_IN_8_10_11_12_13; in.GLYCOGEN_IN1_8_10_11_12_13; in.GLYCOGEN_IN2_8_10_11_12_13; cauchy(x9(6,:),x1(1,:)); cauchy(x1(6,:),x9(5,:)); cauchy(x1(10,:),x9(1,:)); cauchy(x1(7,:),x9(5,:)); cauchy(x1(11,:),x9(1,:))];
x14=A14\(B14*sparse(y14));

%%
A15=zeros(5);
A15(1,1)=-b(25)-f(26)-f(27);
A15(1,2)=f(25);
A15(2,1)=b(25);
A15(2,2)=-b(23)-f(24)-f(25);
A15(2,3)=b(24)*(1-f(44));
A15(3,2)=f(24)*(1-f(44));
A15(3,3)=-b(6)-b(24)-f(27)-f(38);
A15(3,4)=b(27);
A15(4,3)=f(27);
A15(4,4)=-b(27)-b(30)-f(28);
A15(4,5)=f(30);
A15(5,4)=b(30);
A15(5,5)=-b(29)-f(30);

B15=zeros(5,8);
B15(1,1)=-b(26);
B15(1,2)=-b(27);
B15(2,3)=-f(23);
B15(2,4)=-b(24)*f(44);
B15(3,5)=-f(24)*f(44);
B15(3,6)=-f(6);
B15(4,7)=-b(28);
B15(5,8)=-f(29);

y15=[ cauchy(x2(4,:),x5(5,:)); cauchy(x2(5,:),x5(5,:)); cauchy(x1(1,:),x11(2,:)); cauchy(x11(3,:),x1(1,:)); cauchy(x11(4,:),x1(1,:)); in.RIBOSE_IN_6_7_9_10_11; cauchy(x1(5,:),x12(1,:)); cauchy(x1(7,:),x12(1,:))];
x15=A15\(B15*sparse(y15));

%%
A16=zeros(4);
A16(1,1)=-b(20)-f(21)-f(23);
A16(1,2)=f(20);
A16(2,1)=b(20);
A16(2,2)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A16(2,3)=b(10)*(1-f(43));
A16(3,2)=f(10)*(1-f(43));
A16(3,3)=-b(10)-b(26)-b(28)-f(11)-f(36);
A16(3,4)=b(11);
A16(4,3)=f(11);
A16(4,4)=-b(11)-f(12);

B16=zeros(4,10);
B16(1,1)=-b(21);
B16(1,2)=-b(23);
B16(2,3)=-b(10)*f(43);
B16(2,4)=-f(7);
B16(2,5)=-f(8);
B16(2,6)=-f(9);
B16(3,7)=-f(10)*f(43);
B16(3,8)=-f(26);
B16(3,9)=-f(28);
B16(4,10)=-b(12);

y16=[ cauchy(x1(1,:),x13(2,:)); cauchy(x1(2,:),x14(1,:)); cauchy(x13(4,:),x1(1,:)); in.GLC_IN_8_9_10_11_12_13; in.GLYCOGEN_IN1_8_9_10_11_12_13; in.GLYCOGEN_IN2_8_9_10_11_12_13; cauchy(x13(5,:),x1(1,:)); cauchy(x1(6,:),x13(1,:)); cauchy(x4(1,:),x9(1,:)); cauchy(x4(2,:),x9(1,:))];
x16=A16\(B16*sparse(y16));

%%
A17=zeros(1);
A17(1,1)=-b(21)-f(22);

B17=zeros(1,2);
B17(1,1)=-f(21);
B17(1,2)=-b(22);

y17=[ cauchy(x1(1,:),x13(3,:)); cauchy(x2(3,:),x9(1,:))];
x17=A17\(B17*sparse(y17));

%%
A18=zeros(4);
A18(1,1)=-b(23)-f(24)-f(25);
A18(1,2)=b(24)*(1-f(44));
A18(2,1)=f(24)*(1-f(44));
A18(2,2)=-b(6)-b(24)-f(27)-f(38);
A18(2,3)=b(27);
A18(3,2)=f(27);
A18(3,3)=-b(27)-b(30)-f(28);
A18(3,4)=f(30);
A18(4,3)=b(30);
A18(4,4)=-b(29)-f(30);

B18=zeros(4,7);
B18(1,1)=-f(23);
B18(1,2)=-b(24)*f(44);
B18(1,3)=-b(25);
B18(2,4)=-f(24)*f(44);
B18(2,5)=-f(6);
B18(3,6)=-b(28);
B18(4,7)=-f(29);

y18=[ cauchy(x1(1,:),x14(2,:)); cauchy(x14(3,:),x1(1,:)); cauchy(x1(1,:),x15(1,:)); cauchy(x14(1,:),x1(1,:)); in.RIBOSE_IN_6_7_8_9_10_11; cauchy(x1(5,:),x13(1,:)); cauchy(x1(7,:),x13(1,:))];
x18=A18\(B18*sparse(y18));

%%
A19=zeros(1);
A19(1,1)=-b(25)-f(26)-f(27);

B19=zeros(1,3);
B19(1,1)=-f(25);
B19(1,2)=-b(26);
B19(1,3)=-b(27);

y19=[ cauchy(x1(1,:),x15(2,:)); cauchy(x2(4,:),x9(1,:)); cauchy(x2(5,:),x9(1,:))];
x19=A19\(B19*sparse(y19));

%%
A20=zeros(3);
A20(1,1)=-b(10)-b(26)-b(28)-f(11)-f(36);
A20(1,2)=f(10);
A20(1,3)=b(11);
A20(2,1)=b(10);
A20(2,2)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A20(3,1)=f(11);
A20(3,3)=-b(11)-f(12);

B20=zeros(3,7);
B20(1,1)=-f(26);
B20(1,2)=-f(28);
B20(2,3)=-b(20);
B20(2,5)=-f(7);
B20(2,6)=-f(8);
B20(2,7)=-f(9);
B20(3,4)=-b(12);

y20=[ cauchy(x1(4,:),x13(1,:)); cauchy(x3(1,:),x9(1,:)); cauchy(x1(2,:),x13(3,:)); cauchy(x3(2,:),x9(1,:)); in.GLC_IN_7_9_10_11_12_13; in.GLYCOGEN_IN1_7_9_10_11_12_13; in.GLYCOGEN_IN2_7_9_10_11_12_13];
x20=A20\(B20*sparse(y20));

%%
A21=zeros(3);
A21(1,1)=-b(7)-b(8)-b(9)-f(10)-f(20)-f(35);
A21(1,2)=b(10)*(1-f(43));
A21(2,1)=f(10)*(1-f(43));
A21(2,2)=-b(10)-b(26)-b(28)-f(11)-f(36);
A21(2,3)=b(11);
A21(3,2)=f(11);
A21(3,3)=-b(11)-f(12);

B21=zeros(3,9);
B21(1,1)=-b(10)*f(43);
B21(1,2)=-b(20);
B21(1,7)=-f(7);
B21(1,8)=-f(8);
B21(1,9)=-f(9);
B21(2,3)=-f(10)*f(43);
B21(2,4)=-f(26);
B21(2,5)=-f(28);
B21(3,6)=-b(12);

y21=[ cauchy(x20(1,:),x1(1,:)); cauchy(x1(2,:),x16(1,:)); cauchy(x20(2,:),x1(1,:)); cauchy(x2(2,:),x13(1,:)); cauchy(x7(3,:),x9(1,:)); cauchy(x7(2,:),x9(1,:)); in.GLC_IN_7_8_9_10_11_12_13; in.GLYCOGEN_IN1_7_8_9_10_11_12_13; in.GLYCOGEN_IN2_7_8_9_10_11_12_13];
x21=A21\(B21*sparse(y21));

%%
A22=zeros(2);
A22(1,1)=-b(27)-b(30)-f(28);
A22(1,2)=f(30);
A22(2,1)=b(30);
A22(2,2)=-b(29)-f(30);

B22=zeros(2,3);
B22(1,1)=-f(27);
B22(1,2)=-b(28);
B22(2,3)=-f(29);

y22=[ cauchy(x2(2,:),x18(2,:)); cauchy(x7(1,:),x13(1,:)); cauchy(x7(2,:),x13(1,:))];
x22=A22\(B22*sparse(y22));

%%
mol.g6p=x21(1,:);
mol.f6p=x21(2,:);
mol.fbp=x21(3,:);
mol.gap=x9(1,:);
mol.dhap=x10(1,:);
mol.bpg13=x5(1,:);
mol.pg3=x5(2,:);
mol.bpg23=x5(3,:);
mol.pg2=x5(4,:);
mol.pep=x2(1,:);
mol.pyr=x6(1,:);
mol.pg6=x16(1,:);
mol.ddg6p=x17(1,:);
mol.ru5p=x18(1,:);
mol.r5p=x18(2,:);
mol.xu5p=x19(1,:);
mol.s7p=x22(1,:);
mol.sbp=x22(2,:);
mol.e4p=x13(1,:);
mol.H=x1(1,:);
mol.NADPH=x1(2,:);
mol.NADH=x1(3,:);