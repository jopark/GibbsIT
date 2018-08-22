function mol=ccell_H(free_net,free_xch,in)

kernel_net=[ -0 -3 3 1 -0 -0 -0 -0 -0 -0 -1 1 -0 -0;
 -0 2.5 -2.5 -0 1 -0 -0 -0 -0 0.5 1 -0.5 0.5 -0.5;
 -0 -6 6 -0 -0 1 -0 -0 -0 -0 -2 2 -0 -0;
 -0 -0 -0 -0 -0 -0 1 -0 -0 -0 -0 -0 -1 -1;
 -1 0.5 -0.5 -0 -0 -0 -0 1 1 0.5 1 0.5 0.5 0.5;
 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
 -0 -2.5 2.5 -0 -0 -0 -0 -0 1 0.5 -0 1.5 0.5 0.5;
 -0 -0.5 1.5 -0 -0 -0 -0 -0 -0 0.5 -0 0.5 0.5 0.5;
 -0 -0.5 1.5 -0 -0 -0 -0 -0 -0 0.5 -0 0.5 0.5 0.5;
 -0 -0.5 0.5 -0 -0 -0 -0 -0 -0 -0.5 -0 0.5 0.5 0.5;
 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1 1;
 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1 1;
 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1;
 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1;
 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 1;
 -0 3 -3 -0 -0 -0 -0 -0 -0 -0 1 -1 -0 -0;
 -0 3 -3 -0 -0 -0 -0 -0 -0 -0 1 -1 -0 -0;
 -0 1 -1 -0 -0 -0 -0 -0 -0 -0 1 -0 -0 -0;
 -0 2 -2 -0 -0 -0 -0 -0 -0 -0 -0 -1 -0 -0;
 -0 1 -1 0 0 0 0 0 0 0 0 -1 0 0;
 -0 1 -1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0;
 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
 -0 -0 1 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0 -0;
 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 1;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

kernel_xch=[ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0;
 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];

net=sparse(kernel_net)*free_net;
xch=sparse(kernel_xch)*free_xch;

f=xch+max(0,net);
b=xch-min(0,net);

%%
A1=zeros(91);
A1(1,1)=-b(2)-b(10)-b(15)-b(19)-b(23)-f(9)-f(10)-f(17)-f(19)-f(26)-f(7)*f(36)-b(18)*f(37)-b(7)*f(36)-f(18)*f(37);
A1(1,7)=f(10);
A1(1,13)=b(7)*f(36);
A1(1,23)=f(15);
A1(1,24)=f(19);
A1(1,25)=f(23)+b(9);
A1(1,26)=b(10);
A1(1,27)=b(17)+f(18)*f(37);
A1(1,28)=b(19);
A1(1,36)=f(7)*f(36);
A1(1,55)=b(18)*f(37);
A1(2,2)=-b(3)-b(16)-b(17)-f(27);
A1(2,16)=f(16);
A1(2,19)=f(17);
A1(3,3)=-b(4)-b(11)-f(28);
A1(3,29)=f(11);
A1(4,4)=-b(19)-f(20)-f(21);
A1(4,11)=b(20);
A1(4,12)=b(21);
A1(4,20)=f(19);
A1(5,5)=-b(7)-b(20)-b(22)-f(8)-f(30);
A1(5,8)=f(22);
A1(5,14)=f(20);
A1(5,22)=f(7);
A1(5,30)=b(8);
A1(6,6)=-b(19)-f(20)-f(21);
A1(6,9)=b(21);
A1(6,13)=b(20);
A1(6,27)=f(19);
A1(7,1)=b(10);
A1(7,7)=-b(9)-f(10)-f(23)-f(31);
A1(7,30)=f(9);
A1(7,31)=b(23);
A1(8,5)=b(22);
A1(8,8)=-b(21)-b(24)-f(22);
A1(8,15)=f(21);
A1(8,31)=f(24);
A1(9,6)=f(21);
A1(9,9)=-b(21)-b(24)-f(22);
A1(9,13)=b(22);
A1(9,32)=f(24);
A1(10,10)=-b(9)-f(10)-f(23)-f(31);
A1(10,18)=b(10);
A1(10,32)=b(23);
A1(10,33)=f(9);
A1(11,4)=f(20);
A1(11,11)=-b(7)-b(20)-b(22)-f(8)-f(30);
A1(11,12)=f(22);
A1(11,16)=f(7);
A1(11,34)=b(8);
A1(12,4)=f(21);
A1(12,11)=b(22);
A1(12,12)=-b(21)-b(24)-f(22);
A1(12,35)=f(24);
A1(13,1)=f(7)*f(36);
A1(13,6)=f(20);
A1(13,9)=f(22);
A1(13,13)=-b(7)-b(20)-b(22)-f(8)-f(30);
A1(13,33)=b(8);
A1(13,36)=f(7)*(1-f(36));
A1(14,5)=b(20);
A1(14,14)=-b(22)-f(20)-f(23)-f(33);
A1(14,37)=f(22);
A1(14,38)=b(23);
A1(15,8)=b(21);
A1(15,15)=-b(18)-f(21)-f(32);
A1(15,20)=f(18);
A1(16,2)=b(16);
A1(16,11)=b(7);
A1(16,16)=-b(5)-b(6)-f(7)-f(16)-f(29);
A1(17,17)=-b(9)-b(10)-b(20)-b(21)-f(11)-f(22);
A1(17,39)=f(9);
A1(17,40)=f(10);
A1(17,41)=f(20)+f(21);
A1(17,42)=b(11);
A1(17,43)=b(22);
A1(18,10)=f(10);
A1(18,18)=-b(9)-b(10)-b(20)-b(21)-f(11)-f(22);
A1(18,44)=f(9);
A1(18,45)=f(20)+f(21);
A1(18,46)=b(11);
A1(18,47)=b(22);
A1(19,2)=b(17);
A1(19,19)=-b(16)-f(17);
A1(19,22)=f(16);
A1(20,4)=b(19);
A1(20,15)=b(18);
A1(20,20)=-b(17)-f(18)-f(19);
A1(20,21)=f(17);
A1(21,20)=b(17);
A1(21,21)=-b(16)-f(17);
A1(21,36)=f(16);
A1(22,5)=b(7);
A1(22,19)=b(16);
A1(22,22)=-b(5)-b(6)-f(7)-f(16)-f(29);
A1(23,1)=b(15);
A1(23,23)=-b(14)-f(15);
A1(23,48)=f(14);
A1(24,1)=b(19);
A1(24,24)=-b(17)-f(18)-f(19);
A1(24,49)=f(17);
A1(24,50)=b(18);
A1(25,1)=f(9)+b(23);
A1(25,25)=-b(9)-f(10)-f(23)-f(31);
A1(25,29)=b(10);
A1(26,1)=f(10);
A1(26,26)=-b(9)-b(10)-b(20)-b(21)-f(11)-f(22);
A1(26,51)=f(9);
A1(26,52)=f(20)+f(21);
A1(26,53)=b(11);
A1(26,54)=b(22);
A1(27,1)=f(17)+b(18)*f(37);
A1(27,6)=b(19);
A1(27,27)=-b(17)-f(18)-f(19);
A1(27,55)=b(18)*(1-f(37));
A1(28,1)=f(19);
A1(28,28)=-b(19)-f(20)-f(21);
A1(28,29)=b(20)+b(21);
A1(29,3)=b(11);
A1(29,25)=f(10);
A1(29,28)=f(20)+f(21);
A1(29,29)=-b(9)-b(10)-b(20)-b(21)-f(11)-f(22);
A1(29,56)=f(9);
A1(29,57)=b(22);
A1(30,5)=f(8);
A1(30,7)=b(9);
A1(30,30)=-b(8)-f(9);
A1(31,7)=f(23);
A1(31,8)=b(24);
A1(31,31)=-b(23)-f(24);
A1(32,9)=b(24);
A1(32,10)=f(23);
A1(32,32)=-b(23)-f(24);
A1(33,10)=b(9);
A1(33,13)=f(8);
A1(33,33)=-b(8)-f(9);
A1(34,11)=f(8);
A1(34,34)=-b(8)-f(9);
A1(34,40)=b(9);
A1(35,12)=b(24);
A1(35,35)=-b(23)-f(24);
A1(35,40)=f(23);
A1(36,1)=b(7)*f(36);
A1(36,13)=b(7)*(1-f(36));
A1(36,21)=b(16);
A1(36,36)=-b(5)-b(6)-f(7)-f(16)-f(29);
A1(37,14)=b(22);
A1(37,37)=-b(21)-b(24)-f(22);
A1(37,38)=f(24);
A1(37,55)=f(21);
A1(38,14)=f(23);
A1(38,37)=b(24);
A1(38,38)=-b(23)-f(24);
A1(39,17)=b(9);
A1(39,39)=-b(8)-f(9);
A1(39,43)=f(8);
A1(40,17)=b(10);
A1(40,34)=f(9);
A1(40,35)=b(23);
A1(40,40)=-b(9)-f(10)-f(23)-f(31);
A1(41,17)=b(20)+b(21);
A1(41,41)=-b(19)-f(20)-f(21);
A1(41,58)=f(19);
A1(42,17)=f(11);
A1(42,42)=-b(11)-f(12);
A1(42,59)=b(12);
A1(43,17)=f(22);
A1(43,39)=b(8);
A1(43,43)=-b(7)-b(20)-b(22)-f(8)-f(30);
A1(43,60)=f(7);
A1(43,61)=f(20);
A1(44,18)=b(9);
A1(44,44)=-b(8)-f(9);
A1(44,47)=f(8);
A1(45,18)=b(20)+b(21);
A1(45,45)=-b(19)-f(20)-f(21);
A1(45,62)=f(19);
A1(46,18)=f(11);
A1(46,46)=-b(11)-f(12);
A1(46,63)=b(12);
A1(47,18)=f(22);
A1(47,44)=b(8);
A1(47,47)=-b(7)-b(20)-b(22)-f(8)-f(30);
A1(47,64)=f(7);
A1(47,65)=f(20);
A1(48,23)=b(14);
A1(48,48)=-b(13)-f(14);
A1(48,66)=f(13);
A1(49,24)=b(17);
A1(49,49)=-b(16)-f(17);
A1(49,67)=f(16);
A1(50,24)=f(18);
A1(50,50)=-b(18)-f(21)-f(32);
A1(50,68)=b(21);
A1(51,26)=b(9);
A1(51,51)=-b(8)-f(9);
A1(51,54)=f(8);
A1(52,26)=b(20)+b(21);
A1(52,52)=-b(19)-f(20)-f(21);
A1(52,69)=f(19);
A1(53,26)=f(11);
A1(53,53)=-b(11)-f(12);
A1(53,66)=b(12);
A1(54,26)=f(22);
A1(54,51)=b(8);
A1(54,54)=-b(7)-b(20)-b(22)-f(8)-f(30);
A1(54,70)=f(7);
A1(54,71)=f(20);
A1(55,1)=f(18)*f(37);
A1(55,27)=f(18)*(1-f(37));
A1(55,37)=b(21);
A1(55,55)=-b(18)-f(21)-f(32);
A1(56,29)=b(9);
A1(56,56)=-b(8)-f(9);
A1(56,57)=f(8);
A1(57,29)=f(22);
A1(57,56)=b(8);
A1(57,57)=-b(7)-b(20)-b(22)-f(8)-f(30);
A1(57,67)=f(7);
A1(57,72)=f(20);
A1(58,41)=b(19);
A1(58,58)=-b(17)-f(18)-f(19);
A1(58,73)=f(17);
A1(58,74)=b(18);
A1(59,42)=f(12);
A1(59,59)=-b(12)-f(13)-f(34);
A1(59,75)=b(13);
A1(60,43)=b(7);
A1(60,60)=-b(5)-b(6)-f(7)-f(16)-f(29);
A1(60,73)=b(16);
A1(61,43)=b(20);
A1(61,61)=-b(22)-f(20)-f(23)-f(33);
A1(61,76)=f(22);
A1(61,77)=b(23);
A1(62,45)=b(19);
A1(62,62)=-b(17)-f(18)-f(19);
A1(62,78)=f(17);
A1(62,79)=b(18);
A1(63,46)=f(12);
A1(63,63)=-b(12)-f(13)-f(34);
A1(63,80)=b(13);
A1(64,47)=b(7);
A1(64,64)=-b(5)-b(6)-f(7)-f(16)-f(29);
A1(64,78)=b(16);
A1(65,47)=b(20);
A1(65,65)=-b(22)-f(20)-f(23)-f(33);
A1(65,81)=f(22);
A1(65,82)=b(23);
A1(66,48)=b(13);
A1(66,53)=f(12);
A1(66,66)=-b(12)-f(13)-f(34);
A1(67,49)=b(16);
A1(67,57)=b(7);
A1(67,67)=-b(5)-b(6)-f(7)-f(16)-f(29);
A1(68,50)=f(21);
A1(68,68)=-b(21)-b(24)-f(22);
A1(68,72)=b(22);
A1(68,83)=f(24);
A1(69,52)=b(19);
A1(69,69)=-b(17)-f(18)-f(19);
A1(69,84)=f(17);
A1(69,85)=b(18);
A1(70,54)=b(7);
A1(70,70)=-b(5)-b(6)-f(7)-f(16)-f(29);
A1(70,84)=b(16);
A1(71,54)=b(20);
A1(71,71)=-b(22)-f(20)-f(23)-f(33);
A1(71,86)=f(22);
A1(71,87)=b(23);
A1(72,57)=b(20);
A1(72,68)=f(22);
A1(72,72)=-b(22)-f(20)-f(23)-f(33);
A1(72,83)=b(23);
A1(73,58)=b(17);
A1(73,60)=f(16);
A1(73,73)=-b(16)-f(17);
A1(74,58)=f(18);
A1(74,74)=-b(18)-f(21)-f(32);
A1(74,76)=b(21);
A1(75,59)=f(13);
A1(75,75)=-b(13)-f(14);
A1(75,88)=b(14);
A1(76,61)=b(22);
A1(76,74)=f(21);
A1(76,76)=-b(21)-b(24)-f(22);
A1(76,77)=f(24);
A1(77,61)=f(23);
A1(77,76)=b(24);
A1(77,77)=-b(23)-f(24);
A1(78,62)=b(17);
A1(78,64)=f(16);
A1(78,78)=-b(16)-f(17);
A1(79,62)=f(18);
A1(79,79)=-b(18)-f(21)-f(32);
A1(79,81)=b(21);
A1(80,63)=f(13);
A1(80,80)=-b(13)-f(14);
A1(80,89)=b(14);
A1(81,65)=b(22);
A1(81,79)=f(21);
A1(81,81)=-b(21)-b(24)-f(22);
A1(81,82)=f(24);
A1(82,65)=f(23);
A1(82,81)=b(24);
A1(82,82)=-b(23)-f(24);
A1(83,68)=b(24);
A1(83,72)=f(23);
A1(83,83)=-b(23)-f(24);
A1(84,69)=b(17);
A1(84,70)=f(16);
A1(84,84)=-b(16)-f(17);
A1(85,69)=f(18);
A1(85,85)=-b(18)-f(21)-f(32);
A1(85,86)=b(21);
A1(86,71)=b(22);
A1(86,85)=f(21);
A1(86,86)=-b(21)-b(24)-f(22);
A1(86,87)=f(24);
A1(87,71)=f(23);
A1(87,86)=b(24);
A1(87,87)=-b(23)-f(24);
A1(88,75)=f(14);
A1(88,88)=-b(14)-f(15);
A1(88,90)=b(15);
A1(89,80)=f(14);
A1(89,89)=-b(14)-f(15);
A1(89,91)=b(15);
A1(90,88)=f(15);
A1(90,90)=-b(15)-f(35);
A1(91,89)=f(15);
A1(91,91)=-b(15)-f(35);

B1=zeros(91,17);
B1(1,1)=-f(2);
B1(2,2)=-f(3);
B1(3,3)=-f(4);
B1(16,4)=-f(5);
B1(16,5)=-f(6);
B1(22,6)=-f(5);
B1(22,7)=-f(6);
B1(36,8)=-f(5);
B1(36,9)=-f(6);
B1(60,10)=-f(5);
B1(60,11)=-f(6);
B1(64,12)=-f(5);
B1(64,13)=-f(6);
B1(67,14)=-f(5);
B1(67,15)=-f(6);
B1(70,16)=-f(5);
B1(70,17)=-f(6);

y1=[ in.H_IN_1; in.NADPH_IN_1; in.NADH_IN_1; in.GLC_IN_7; in.GLYCOGEN_IN1_7; in.GLC_IN_9; in.GLYCOGEN_IN1_9; in.GLC_IN_8; in.GLYCOGEN_IN1_8; in.GLC_IN_12; in.GLYCOGEN_IN1_12; in.GLC_IN_13; in.GLYCOGEN_IN1_13; in.GLC_IN_10; in.GLYCOGEN_IN1_10; in.GLC_IN_11; in.GLYCOGEN_IN1_11];
x1=sparse(A1)\(B1*sparse(y1));

%%
A2=zeros(27);
A2(1,1)=-b(15)-f(35);
A2(1,12)=f(15);
A2(2,2)=-b(19)-f(20)-f(21);
A2(2,3)=b(20);
A2(2,4)=b(21);
A2(2,13)=f(19);
A2(3,2)=f(20);
A2(3,3)=-b(7)-b(20)-b(22)-f(8)-f(30);
A2(3,4)=f(22);
A2(3,7)=b(8);
A2(3,14)=f(7)*(1-f(36));
A2(4,2)=f(21);
A2(4,3)=b(22);
A2(4,4)=-b(21)-b(24)-f(22);
A2(4,8)=f(24);
A2(5,5)=-b(9)-b(10)-b(20)-b(21)-f(11)-f(22);
A2(5,6)=f(10);
A2(5,9)=b(11);
A2(5,11)=f(20)+f(21);
A2(5,15)=f(9);
A2(5,16)=b(22);
A2(6,5)=b(10);
A2(6,6)=-b(9)-f(10)-f(23)-f(31);
A2(6,7)=f(9);
A2(6,8)=b(23);
A2(7,3)=f(8);
A2(7,6)=b(9);
A2(7,7)=-b(8)-f(9);
A2(8,4)=b(24);
A2(8,6)=f(23);
A2(8,8)=-b(23)-f(24);
A2(9,5)=f(11);
A2(9,9)=-b(11)-f(12);
A2(9,17)=b(12);
A2(10,10)=-b(17)-f(18)-f(19);
A2(10,11)=b(19);
A2(10,18)=f(17);
A2(10,19)=b(18);
A2(11,5)=b(20)+b(21);
A2(11,10)=f(19);
A2(11,11)=-b(19)-f(20)-f(21);
A2(12,1)=b(15);
A2(12,12)=-b(14)-f(15);
A2(12,20)=f(14);
A2(13,2)=b(19);
A2(13,13)=-b(17)-f(18)-f(19);
A2(13,21)=b(18)*(1-f(37));
A2(14,3)=b(7)*(1-f(36));
A2(14,14)=-b(5)-b(6)-f(7)-f(16)-f(29);
A2(15,5)=b(9);
A2(15,15)=-b(8)-f(9);
A2(15,16)=f(8);
A2(16,5)=f(22);
A2(16,15)=b(8);
A2(16,16)=-b(7)-b(20)-b(22)-f(8)-f(30);
A2(16,22)=f(7);
A2(16,23)=f(20);
A2(17,9)=f(12);
A2(17,17)=-b(12)-f(13)-f(34);
A2(17,20)=b(13);
A2(18,10)=b(17);
A2(18,18)=-b(16)-f(17);
A2(18,22)=f(16);
A2(19,10)=f(18);
A2(19,19)=-b(18)-f(21)-f(32);
A2(19,24)=b(21);
A2(20,12)=b(14);
A2(20,17)=f(13);
A2(20,20)=-b(13)-f(14);
A2(21,13)=f(18)*(1-f(37));
A2(21,21)=-b(18)-f(21)-f(32);
A2(21,25)=b(21);
A2(22,16)=b(7);
A2(22,18)=b(16);
A2(22,22)=-b(5)-b(6)-f(7)-f(16)-f(29);
A2(23,16)=b(20);
A2(23,23)=-b(22)-f(20)-f(23)-f(33);
A2(23,24)=f(22);
A2(23,26)=b(23);
A2(24,19)=f(21);
A2(24,23)=b(22);
A2(24,24)=-b(21)-b(24)-f(22);
A2(24,26)=f(24);
A2(25,21)=f(21);
A2(25,25)=-b(21)-b(24)-f(22);
A2(25,27)=f(24);
A2(26,23)=f(23);
A2(26,24)=b(24);
A2(26,26)=-b(23)-f(24);
A2(27,25)=b(24);
A2(27,27)=-b(23)-f(24);

B2=zeros(27,12);
B2(3,1)=-f(7)*f(36);
B2(13,2)=-f(17);
B2(13,3)=-b(18)*f(37);
B2(14,4)=-b(7)*f(36);
B2(14,5)=-b(16);
B2(14,7)=-f(5);
B2(14,8)=-f(6);
B2(21,6)=-f(18)*f(37);
B2(22,10)=-f(5);
B2(22,11)=-f(6);
B2(25,9)=-b(22);
B2(27,12)=-f(23);

y2=[ cauchy(x1(16,:),x1(1,:)); cauchy(x1(1,:),x1(21,:)); cauchy(x1(15,:),x1(1,:)); cauchy(x1(11,:),x1(1,:)); cauchy(x1(2,:),x1(21,:)); cauchy(x1(20,:),x1(1,:)); in.GLC_IN_7_8; in.GLYCOGEN_IN1_7_8; cauchy(x1(5,:),x1(14,:)); in.GLC_IN_12_13; in.GLYCOGEN_IN1_12_13; cauchy(x1(7,:),x1(14,:))];
x2=sparse(A2)\(B2*sparse(y2));

%%
A3=zeros(6);
A3(1,1)=-b(21)-b(24)-f(22);
A3(1,4)=b(22);
A3(1,5)=f(24);
A3(2,2)=-b(9)-f(10)-f(23)-f(31);
A3(2,5)=b(23);
A3(2,6)=f(9);
A3(3,3)=-b(5)-b(6)-f(7)-f(16)-f(29);
A3(3,4)=b(7);
A3(4,1)=f(22);
A3(4,3)=f(7);
A3(4,4)=-b(7)-b(20)-b(22)-f(8)-f(30);
A3(4,6)=b(8);
A3(5,1)=b(24);
A3(5,2)=f(23);
A3(5,5)=-b(23)-f(24);
A3(6,2)=b(9);
A3(6,4)=f(8);
A3(6,6)=-b(8)-f(9);

B3=zeros(6,6);
B3(1,1)=-f(21);
B3(2,2)=-b(10);
B3(3,3)=-b(16);
B3(3,5)=-f(5);
B3(3,6)=-f(6);
B3(4,4)=-f(20);

y3=[ cauchy(x1(4,:),x1(15,:)); cauchy(x1(1,:),x1(17,:)); cauchy(x1(2,:),x1(19,:)); cauchy(x1(4,:),x1(14,:)); in.GLC_IN_7_9; in.GLYCOGEN_IN1_7_9];
x3=A3\(B3*sparse(y3));

%%
A4=zeros(7);
A4(1,1)=-b(21)-b(24)-f(22);
A4(1,4)=f(24);
A4(1,5)=b(22);
A4(2,2)=-b(9)-f(10)-f(23)-f(31);
A4(2,4)=b(23);
A4(2,6)=f(9);
A4(3,3)=-b(16)-f(17);
A4(3,7)=f(16);
A4(4,1)=b(24);
A4(4,2)=f(23);
A4(4,4)=-b(23)-f(24);
A4(5,1)=f(22);
A4(5,5)=-b(7)-b(20)-b(22)-f(8)-f(30);
A4(5,6)=b(8);
A4(5,7)=f(7)*(1-f(36));
A4(6,2)=b(9);
A4(6,5)=f(8);
A4(6,6)=-b(8)-f(9);
A4(7,3)=b(16);
A4(7,5)=b(7)*(1-f(36));
A4(7,7)=-b(5)-b(6)-f(7)-f(16)-f(29);

B4=zeros(7,8);
B4(1,1)=-f(21);
B4(2,2)=-b(10);
B4(3,3)=-b(17);
B4(5,4)=-f(7)*f(36);
B4(5,5)=-f(20);
B4(7,6)=-b(7)*f(36);
B4(7,7)=-f(5);
B4(7,8)=-f(6);

y4=[ cauchy(x1(6,:),x1(15,:)); cauchy(x1(1,:),x1(18,:)); cauchy(x1(2,:),x1(20,:)); cauchy(x1(22,:),x1(1,:)); cauchy(x1(6,:),x1(14,:)); cauchy(x1(5,:),x1(1,:)); in.GLC_IN_8_9; in.GLYCOGEN_IN1_8_9];
x4=A4\(B4*sparse(y4));

%%
A5=zeros(15);
A5(1,1)=-b(11)-f(12);
A5(1,2)=b(12);
A5(1,5)=f(11);
A5(2,1)=f(12);
A5(2,2)=-b(12)-f(13)-f(34);
A5(2,3)=b(13);
A5(3,2)=f(13);
A5(3,3)=-b(13)-f(14);
A5(3,4)=b(14);
A5(4,3)=f(14);
A5(4,4)=-b(14)-f(15);
A5(5,1)=b(11);
A5(5,5)=-b(9)-b(10)-b(20)-b(21)-f(11)-f(22);
A5(5,6)=f(20)+f(21);
A5(5,8)=b(22);
A5(5,13)=f(9);
A5(6,5)=b(20)+b(21);
A5(6,6)=-b(19)-f(20)-f(21);
A5(6,7)=f(19);
A5(7,6)=b(19);
A5(7,7)=-b(17)-f(18)-f(19);
A5(7,10)=f(17);
A5(7,11)=b(18);
A5(8,5)=f(22);
A5(8,8)=-b(7)-b(20)-b(22)-f(8)-f(30);
A5(8,9)=f(20);
A5(8,12)=f(7);
A5(8,13)=b(8);
A5(9,8)=b(20);
A5(9,9)=-b(22)-f(20)-f(23)-f(33);
A5(9,14)=f(22);
A5(9,15)=b(23);
A5(10,7)=b(17);
A5(10,10)=-b(16)-f(17);
A5(10,12)=f(16);
A5(11,7)=f(18);
A5(11,11)=-b(18)-f(21)-f(32);
A5(11,14)=b(21);
A5(12,8)=b(7);
A5(12,10)=b(16);
A5(12,12)=-b(5)-b(6)-f(7)-f(16)-f(29);
A5(13,5)=b(9);
A5(13,8)=f(8);
A5(13,13)=-b(8)-f(9);
A5(14,9)=b(22);
A5(14,11)=f(21);
A5(14,14)=-b(21)-b(24)-f(22);
A5(14,15)=f(24);
A5(15,9)=f(23);
A5(15,14)=b(24);
A5(15,15)=-b(23)-f(24);

B5=zeros(15,4);
B5(4,1)=-b(15);
B5(5,2)=-f(10);
B5(12,3)=-f(5);
B5(12,4)=-f(6);

y5=[ cauchy(x1(1,:),x2(1,:)); cauchy(x1(1,:),x2(6,:)); in.GLC_IN_11_12_13; in.GLYCOGEN_IN1_11_12_13];
x5=sparse(A5)\(B5*sparse(y5));

%%
A6=zeros(6);
A6(1,1)=-b(7)-b(20)-b(22)-f(8)-f(30);
A6(1,3)=f(22);
A6(1,4)=b(8);
A6(1,6)=f(7)*(1-f(36));
A6(2,2)=-b(9)-f(10)-f(23)-f(31);
A6(2,4)=f(9);
A6(2,5)=b(23);
A6(3,1)=b(22);
A6(3,3)=-b(21)-b(24)-f(22);
A6(3,5)=f(24);
A6(4,1)=f(8);
A6(4,2)=b(9);
A6(4,4)=-b(8)-f(9);
A6(5,2)=f(23);
A6(5,3)=b(24);
A6(5,5)=-b(23)-f(24);
A6(6,1)=b(7)*(1-f(36));
A6(6,6)=-b(5)-b(6)-f(7)-f(16)-f(29);

B6=zeros(6,8);
B6(1,1)=-f(7)*f(36);
B6(1,2)=-f(20);
B6(2,3)=-b(10);
B6(3,4)=-f(21);
B6(6,5)=-b(7)*f(36);
B6(6,6)=-b(16);
B6(6,7)=-f(5);
B6(6,8)=-f(6);

y6=[ cauchy(x3(3,:),x1(1,:)); cauchy(x2(2,:),x1(14,:)); cauchy(x1(1,:),x2(5,:)); cauchy(x2(2,:),x1(15,:)); cauchy(x3(4,:),x1(1,:)); cauchy(x1(2,:),x4(3,:)); in.GLC_IN_7_8_9; in.GLYCOGEN_IN1_7_8_9];
x6=A6\(B6*sparse(y6));

%%
A7=zeros(12);
A7(1,1)=-b(9)-f(10)-f(23)-f(31);
A7(1,2)=b(10);
A7(2,1)=f(10);
A7(2,2)=-b(9)-b(10)-b(20)-b(21)-f(11)-f(22);
A7(2,3)=f(9);
A7(2,4)=f(20)+f(21);
A7(2,5)=b(22);
A7(3,2)=b(9);
A7(3,3)=-b(8)-f(9);
A7(3,5)=f(8);
A7(4,2)=b(20)+b(21);
A7(4,4)=-b(19)-f(20)-f(21);
A7(5,2)=f(22);
A7(5,3)=b(8);
A7(5,5)=-b(7)-b(20)-b(22)-f(8)-f(30);
A7(5,6)=f(7);
A7(5,7)=f(20);
A7(6,5)=b(7);
A7(6,6)=-b(5)-b(6)-f(7)-f(16)-f(29);
A7(6,8)=b(16);
A7(7,5)=b(20);
A7(7,7)=-b(22)-f(20)-f(23)-f(33);
A7(7,9)=f(22);
A7(7,10)=b(23);
A7(8,6)=f(16);
A7(8,8)=-b(16)-f(17);
A7(8,11)=b(17);
A7(9,7)=b(22);
A7(9,9)=-b(21)-b(24)-f(22);
A7(9,10)=f(24);
A7(9,12)=f(21);
A7(10,7)=f(23);
A7(10,9)=b(24);
A7(10,10)=-b(23)-f(24);
A7(11,8)=f(17);
A7(11,11)=-b(17)-f(18)-f(19);
A7(11,12)=b(18);
A7(12,9)=b(21);
A7(12,11)=f(18);
A7(12,12)=-b(18)-f(21)-f(32);

B7=zeros(12,7);
B7(1,1)=-f(9);
B7(1,2)=-b(23);
B7(2,3)=-b(11);
B7(4,4)=-f(19);
B7(6,5)=-f(5);
B7(6,6)=-f(6);
B7(11,7)=-b(19);

y7=[ cauchy(x1(1,:),x2(7,:)); cauchy(x1(1,:),x2(8,:)); cauchy(x1(3,:),x2(9,:)); cauchy(x1(1,:),x2(10,:)); in.GLC_IN_10_12_13; in.GLYCOGEN_IN1_10_12_13; cauchy(x1(1,:),x2(11,:))];
x7=sparse(A7)\(B7*sparse(y7));

%%
A8=zeros(11);
A8(1,1)=-b(9)-b(10)-b(20)-b(21)-f(11)-f(22);
A8(1,3)=b(22);
A8(1,8)=f(9);
A8(1,9)=f(20)+f(21);
A8(2,2)=-b(17)-f(18)-f(19);
A8(2,6)=f(17);
A8(2,7)=b(18);
A8(3,1)=f(22);
A8(3,3)=-b(7)-b(20)-b(22)-f(8)-f(30);
A8(3,4)=f(20);
A8(3,5)=f(7);
A8(3,8)=b(8);
A8(4,3)=b(20);
A8(4,4)=-b(22)-f(20)-f(23)-f(33);
A8(4,10)=f(22);
A8(4,11)=b(23);
A8(5,3)=b(7);
A8(5,5)=-b(5)-b(6)-f(7)-f(16)-f(29);
A8(5,6)=b(16);
A8(6,2)=b(17);
A8(6,5)=f(16);
A8(6,6)=-b(16)-f(17);
A8(7,2)=f(18);
A8(7,7)=-b(18)-f(21)-f(32);
A8(7,10)=b(21);
A8(8,1)=b(9);
A8(8,3)=f(8);
A8(8,8)=-b(8)-f(9);
A8(9,1)=b(20)+b(21);
A8(9,9)=-b(19)-f(20)-f(21);
A8(10,4)=b(22);
A8(10,7)=f(21);
A8(10,10)=-b(21)-b(24)-f(22);
A8(10,11)=f(24);
A8(11,4)=f(23);
A8(11,10)=b(24);
A8(11,11)=-b(23)-f(24);

B8=zeros(11,6);
B8(1,1)=-f(10);
B8(1,2)=-b(11);
B8(2,3)=-b(19);
B8(5,5)=-f(5);
B8(5,6)=-f(6);
B8(9,4)=-f(19);

y8=[ cauchy(x1(1,:),x7(1,:)); cauchy(x1(3,:),x5(1,:)); cauchy(x1(1,:),x5(6,:)); cauchy(x1(1,:),x5(7,:)); in.GLC_IN_10_11_12_13; in.GLYCOGEN_IN1_10_11_12_13];
x8=sparse(A8)\(B8*sparse(y8));

%%
A9=zeros(1);
A9(1,1)=-b(9)-f(10)-f(23)-f(31);

B9=zeros(1,3);
B9(1,1)=-f(9);
B9(1,2)=-b(10);
B9(1,3)=-b(23);

y9=[ cauchy(x1(1,:),x6(4,:)); cauchy(x1(1,:),x7(2,:)); cauchy(x1(1,:),x6(5,:))];
x9=A9\(B9*sparse(y9));

%%
A10=zeros(9);
A10(1,1)=-b(19)-f(20)-f(21);
A10(1,4)=f(19);
A10(2,2)=-b(16)-f(17);
A10(2,4)=b(17);
A10(2,5)=f(16);
A10(3,3)=-b(18)-f(21)-f(32);
A10(3,4)=f(18);
A10(3,6)=b(21);
A10(4,1)=b(19);
A10(4,2)=f(17);
A10(4,3)=b(18);
A10(4,4)=-b(17)-f(18)-f(19);
A10(5,2)=b(16);
A10(5,5)=-b(5)-b(6)-f(7)-f(16)-f(29);
A10(5,7)=b(7)*(1-f(36));
A10(6,3)=f(21);
A10(6,6)=-b(21)-b(24)-f(22);
A10(6,8)=f(24);
A10(7,5)=f(7)*(1-f(36));
A10(7,7)=-b(7)-b(20)-b(22)-f(8)-f(30);
A10(7,9)=b(8);
A10(8,6)=b(24);
A10(8,8)=-b(23)-f(24);
A10(9,7)=f(8);
A10(9,9)=-b(8)-f(9);

B10=zeros(9,11);
B10(1,1)=-b(20);
B10(1,2)=-b(21);
B10(5,3)=-b(7)*f(36);
B10(5,5)=-f(5);
B10(5,6)=-f(6);
B10(6,4)=-b(22);
B10(7,7)=-f(7)*f(36);
B10(7,8)=-f(20);
B10(7,9)=-f(22);
B10(8,10)=-f(23);
B10(9,11)=-b(9);

y10=[ cauchy(x1(11,:),x5(5,:)); cauchy(x1(12,:),x5(5,:)); cauchy(x5(8,:),x1(1,:)); cauchy(x1(5,:),x5(9,:)); in.GLC_IN_8_11_12_13; in.GLYCOGEN_IN1_8_11_12_13; cauchy(x5(12,:),x1(1,:)); cauchy(x1(6,:),x5(9,:)); cauchy(x1(9,:),x5(5,:)); cauchy(x1(7,:),x5(9,:)); cauchy(x1(10,:),x5(5,:))];
x10=A10\(B10*sparse(y10));

%%
A11=zeros(10);
A11(1,1)=-b(22)-f(20)-f(23)-f(33);
A11(1,3)=f(22);
A11(1,4)=b(20);
A11(1,5)=b(23);
A11(2,2)=-b(19)-f(20)-f(21);
A11(2,6)=f(19);
A11(3,1)=b(22);
A11(3,3)=-b(21)-b(24)-f(22);
A11(3,5)=f(24);
A11(3,7)=f(21);
A11(4,1)=f(20);
A11(4,4)=-b(7)-b(20)-b(22)-f(8)-f(30);
A11(4,8)=f(7);
A11(4,9)=b(8);
A11(5,1)=f(23);
A11(5,3)=b(24);
A11(5,5)=-b(23)-f(24);
A11(6,2)=b(19);
A11(6,6)=-b(17)-f(18)-f(19);
A11(6,7)=b(18)*(1-f(37));
A11(7,3)=b(21);
A11(7,6)=f(18)*(1-f(37));
A11(7,7)=-b(18)-f(21)-f(32);
A11(8,4)=b(7);
A11(8,8)=-b(5)-b(6)-f(7)-f(16)-f(29);
A11(8,10)=b(16);
A11(9,4)=f(8);
A11(9,9)=-b(8)-f(9);
A11(10,8)=f(16);
A11(10,10)=-b(16)-f(17);

B11=zeros(10,10);
B11(2,1)=-b(20);
B11(2,2)=-b(21);
B11(4,3)=-f(22);
B11(6,4)=-f(17);
B11(6,5)=-b(18)*f(37);
B11(7,6)=-f(18)*f(37);
B11(8,8)=-f(5);
B11(8,9)=-f(6);
B11(9,7)=-b(9);
B11(10,10)=-b(17);

y11=[ cauchy(x1(13,:),x5(5,:)); cauchy(x1(9,:),x5(5,:)); cauchy(x1(8,:),x5(5,:)); cauchy(x1(1,:),x5(10,:)); cauchy(x5(11,:),x1(1,:)); cauchy(x5(7,:),x1(1,:)); cauchy(x1(7,:),x5(5,:)); in.GLC_IN_9_11_12_13; in.GLYCOGEN_IN1_9_11_12_13; cauchy(x1(2,:),x5(7,:))];
x11=A11\(B11*sparse(y11));

%%
A12=zeros(9);
A12(1,1)=-b(22)-f(20)-f(23)-f(33);
A12(1,3)=b(20);
A12(1,5)=f(22);
A12(1,6)=b(23);
A12(2,2)=-b(16)-f(17);
A12(2,4)=f(16);
A12(3,1)=f(20);
A12(3,3)=-b(7)-b(20)-b(22)-f(8)-f(30);
A12(3,4)=f(7);
A12(3,7)=b(8);
A12(4,2)=b(16);
A12(4,3)=b(7);
A12(4,4)=-b(5)-b(6)-f(7)-f(16)-f(29);
A12(5,1)=b(22);
A12(5,5)=-b(21)-b(24)-f(22);
A12(5,6)=f(24);
A12(5,8)=f(21);
A12(6,1)=f(23);
A12(6,5)=b(24);
A12(6,6)=-b(23)-f(24);
A12(7,3)=f(8);
A12(7,7)=-b(8)-f(9);
A12(8,5)=b(21);
A12(8,8)=-b(18)-f(21)-f(32);
A12(8,9)=f(18)*(1-f(37));
A12(9,8)=b(18)*(1-f(37));
A12(9,9)=-b(17)-f(18)-f(19);

B12=zeros(9,9);
B12(2,1)=-b(17);
B12(3,2)=-f(22);
B12(4,4)=-f(5);
B12(4,5)=-f(6);
B12(7,3)=-b(9);
B12(8,6)=-f(18)*f(37);
B12(9,7)=-f(17);
B12(9,8)=-b(18)*f(37);
B12(9,9)=-b(19);

y12=[ cauchy(x1(2,:),x8(2,:)); cauchy(x1(8,:),x8(1,:)); cauchy(x1(7,:),x8(1,:)); in.GLC_IN_9_10_11_12_13; in.GLYCOGEN_IN1_9_10_11_12_13; cauchy(x8(2,:),x1(1,:)); cauchy(x1(1,:),x8(6,:)); cauchy(x8(7,:),x1(1,:)); cauchy(x1(1,:),x11(2,:))];
x12=A12\(B12*sparse(y12));

%%
A13=zeros(8);
A13(1,1)=-b(17)-f(18)-f(19);
A13(1,2)=f(17);
A13(1,3)=b(18);
A13(2,1)=b(17);
A13(2,2)=-b(16)-f(17);
A13(2,4)=f(16);
A13(3,1)=f(18);
A13(3,3)=-b(18)-f(21)-f(32);
A13(3,5)=b(21);
A13(4,2)=b(16);
A13(4,4)=-b(5)-b(6)-f(7)-f(16)-f(29);
A13(4,6)=b(7)*(1-f(36));
A13(5,3)=f(21);
A13(5,5)=-b(21)-b(24)-f(22);
A13(5,7)=f(24);
A13(6,4)=f(7)*(1-f(36));
A13(6,6)=-b(7)-b(20)-b(22)-f(8)-f(30);
A13(6,8)=b(8);
A13(7,5)=b(24);
A13(7,7)=-b(23)-f(24);
A13(8,6)=f(8);
A13(8,8)=-b(8)-f(9);

B13=zeros(8,10);
B13(1,1)=-b(19);
B13(4,2)=-b(7)*f(36);
B13(4,4)=-f(5);
B13(4,5)=-f(6);
B13(5,3)=-b(22);
B13(6,6)=-f(7)*f(36);
B13(6,7)=-f(20);
B13(6,8)=-f(22);
B13(7,9)=-f(23);
B13(8,10)=-b(9);

y13=[ cauchy(x1(1,:),x10(1,:)); cauchy(x8(3,:),x1(1,:)); cauchy(x1(5,:),x8(4,:)); in.GLC_IN_8_10_11_12_13; in.GLYCOGEN_IN1_8_10_11_12_13; cauchy(x8(5,:),x1(1,:)); cauchy(x1(6,:),x8(4,:)); cauchy(x1(9,:),x8(1,:)); cauchy(x1(7,:),x8(4,:)); cauchy(x1(10,:),x8(1,:))];
x13=A13\(B13*sparse(y13));

%%
A14=zeros(5);
A14(1,1)=-b(19)-f(20)-f(21);
A14(1,2)=f(19);
A14(2,1)=b(19);
A14(2,2)=-b(17)-f(18)-f(19);
A14(2,3)=b(18)*(1-f(37));
A14(3,2)=f(18)*(1-f(37));
A14(3,3)=-b(18)-f(21)-f(32);
A14(3,4)=b(21);
A14(4,3)=f(21);
A14(4,4)=-b(21)-b(24)-f(22);
A14(4,5)=f(24);
A14(5,4)=b(24);
A14(5,5)=-b(23)-f(24);

B14=zeros(5,7);
B14(1,1)=-b(20);
B14(1,2)=-b(21);
B14(2,3)=-f(17);
B14(2,4)=-b(18)*f(37);
B14(3,5)=-f(18)*f(37);
B14(4,6)=-b(22);
B14(5,7)=-f(23);

y14=[ cauchy(x2(3,:),x5(5,:)); cauchy(x2(4,:),x5(5,:)); cauchy(x1(1,:),x10(2,:)); cauchy(x10(3,:),x1(1,:)); cauchy(x10(4,:),x1(1,:)); cauchy(x1(5,:),x11(1,:)); cauchy(x1(7,:),x11(1,:))];
x14=A14\(B14*sparse(y14));

%%
A15=zeros(4);
A15(1,1)=-b(16)-f(17);
A15(1,2)=f(16);
A15(2,1)=b(16);
A15(2,2)=-b(5)-b(6)-f(7)-f(16)-f(29);
A15(2,3)=b(7)*(1-f(36));
A15(3,2)=f(7)*(1-f(36));
A15(3,3)=-b(7)-b(20)-b(22)-f(8)-f(30);
A15(3,4)=b(8);
A15(4,3)=f(8);
A15(4,4)=-b(8)-f(9);

B15=zeros(4,8);
B15(1,1)=-b(17);
B15(2,2)=-b(7)*f(36);
B15(2,3)=-f(5);
B15(2,4)=-f(6);
B15(3,5)=-f(7)*f(36);
B15(3,6)=-f(20);
B15(3,7)=-f(22);
B15(4,8)=-b(9);

y15=[ cauchy(x1(2,:),x13(1,:)); cauchy(x12(3,:),x1(1,:)); in.GLC_IN_8_9_10_11_12_13; in.GLYCOGEN_IN1_8_9_10_11_12_13; cauchy(x12(4,:),x1(1,:)); cauchy(x1(6,:),x12(1,:)); cauchy(x4(1,:),x8(1,:)); cauchy(x4(2,:),x8(1,:))];
x15=A15\(B15*sparse(y15));

%%
A16=zeros(4);
A16(1,1)=-b(17)-f(18)-f(19);
A16(1,2)=b(18)*(1-f(37));
A16(2,1)=f(18)*(1-f(37));
A16(2,2)=-b(18)-f(21)-f(32);
A16(2,3)=b(21);
A16(3,2)=f(21);
A16(3,3)=-b(21)-b(24)-f(22);
A16(3,4)=f(24);
A16(4,3)=b(24);
A16(4,4)=-b(23)-f(24);

B16=zeros(4,6);
B16(1,1)=-f(17);
B16(1,2)=-b(18)*f(37);
B16(1,3)=-b(19);
B16(2,4)=-f(18)*f(37);
B16(3,5)=-b(22);
B16(4,6)=-f(23);

y16=[ cauchy(x1(1,:),x13(2,:)); cauchy(x13(3,:),x1(1,:)); cauchy(x1(1,:),x14(1,:)); cauchy(x13(1,:),x1(1,:)); cauchy(x1(5,:),x12(1,:)); cauchy(x1(7,:),x12(1,:))];
x16=A16\(B16*sparse(y16));

%%
A17=zeros(1);
A17(1,1)=-b(19)-f(20)-f(21);

B17=zeros(1,3);
B17(1,1)=-f(19);
B17(1,2)=-b(20);
B17(1,3)=-b(21);

y17=[ cauchy(x1(1,:),x14(2,:)); cauchy(x2(3,:),x8(1,:)); cauchy(x2(4,:),x8(1,:))];
x17=A17\(B17*sparse(y17));

%%
A18=zeros(3);
A18(1,1)=-b(7)-b(20)-b(22)-f(8)-f(30);
A18(1,2)=f(7);
A18(1,3)=b(8);
A18(2,1)=b(7);
A18(2,2)=-b(5)-b(6)-f(7)-f(16)-f(29);
A18(3,1)=f(8);
A18(3,3)=-b(8)-f(9);

B18=zeros(3,6);
B18(1,1)=-f(20);
B18(1,2)=-f(22);
B18(2,3)=-b(16);
B18(2,5)=-f(5);
B18(2,6)=-f(6);
B18(3,4)=-b(9);

y18=[ cauchy(x1(4,:),x12(1,:)); cauchy(x3(1,:),x8(1,:)); cauchy(x1(2,:),x12(2,:)); cauchy(x3(2,:),x8(1,:)); in.GLC_IN_7_9_10_11_12_13; in.GLYCOGEN_IN1_7_9_10_11_12_13];
x18=A18\(B18*sparse(y18));

%%
A19=zeros(3);
A19(1,1)=-b(5)-b(6)-f(7)-f(16)-f(29);
A19(1,2)=b(7)*(1-f(36));
A19(2,1)=f(7)*(1-f(36));
A19(2,2)=-b(7)-b(20)-b(22)-f(8)-f(30);
A19(2,3)=b(8);
A19(3,2)=f(8);
A19(3,3)=-b(8)-f(9);

B19=zeros(3,8);
B19(1,1)=-b(7)*f(36);
B19(1,2)=-b(16);
B19(1,7)=-f(5);
B19(1,8)=-f(6);
B19(2,3)=-f(7)*f(36);
B19(2,4)=-f(20);
B19(2,5)=-f(22);
B19(3,6)=-b(9);

y19=[ cauchy(x18(1,:),x1(1,:)); cauchy(x1(2,:),x15(1,:)); cauchy(x18(2,:),x1(1,:)); cauchy(x2(2,:),x12(1,:)); cauchy(x6(3,:),x8(1,:)); cauchy(x6(2,:),x8(1,:)); in.GLC_IN_7_8_9_10_11_12_13; in.GLYCOGEN_IN1_7_8_9_10_11_12_13];
x19=A19\(B19*sparse(y19));

%%
A20=zeros(2);
A20(1,1)=-b(21)-b(24)-f(22);
A20(1,2)=f(24);
A20(2,1)=b(24);
A20(2,2)=-b(23)-f(24);

B20=zeros(2,3);
B20(1,1)=-f(21);
B20(1,2)=-b(22);
B20(2,3)=-f(23);

y20=[ cauchy(x2(2,:),x16(2,:)); cauchy(x6(1,:),x12(1,:)); cauchy(x6(2,:),x12(1,:))];
x20=A20\(B20*sparse(y20));

%%
mol.g6p=x19(1,:);
mol.f6p=x19(2,:);
mol.fbp=x19(3,:);
mol.gap=x8(1,:);
mol.dhap=x9(1,:);
mol.bpg13=x5(1,:);
mol.pg3=x5(2,:);
mol.bpg23=x5(3,:);
mol.pg2=x5(4,:);
mol.pep=x2(1,:);
mol.pg6=x15(1,:);
mol.ru5p=x16(1,:);
mol.r5p=x16(2,:);
mol.xu5p=x17(1,:);
mol.s7p=x20(1,:);
mol.sbp=x20(2,:);
mol.e4p=x12(1,:);
mol.H=x1(1,:);
mol.NADPH=x1(2,:);
mol.NADH=x1(3,:);