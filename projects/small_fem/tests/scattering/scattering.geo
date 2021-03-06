Include "parameters_gmsh_getdp.dat";

In_n  = Sqrt(eps_re_In);
In_n  = Sqrt(eps_re_In);
Out_n = Sqrt(eps_re_Out);
Out_n = Sqrt(eps_re_Out);

paramaille_pml = paramaille / 1.2;

In_lc         = lambda0 / (paramaille     * In_n);
Out_lc        = lambda0 / (paramaille     * Out_n);
PML_lc        = lambda0 / (paramaille_pml * Out_n);
CenterScat_lc = lambda0 / (paramaille     * In_n);//lambda0/(paramaille*2);

Point(1)  = {-dom_x/2, -dom_y/2, -dom_z/2, Out_lc};
Point(2)  = {-dom_x/2,  dom_y/2, -dom_z/2, Out_lc};
Point(3)  = { dom_x/2,  dom_y/2, -dom_z/2, Out_lc};
Point(4)  = { dom_x/2, -dom_y/2, -dom_z/2, Out_lc};
Point(5)  = {-dom_x/2, -dom_y/2,  dom_z/2, Out_lc};
Point(6)  = {-dom_x/2,  dom_y/2,  dom_z/2, Out_lc};
Point(7)  = { dom_x/2,  dom_y/2,  dom_z/2, Out_lc};
Point(8)  = { dom_x/2, -dom_y/2,  dom_z/2, Out_lc};

// PML top and bot
Point(9)   = {-dom_x/2, -dom_y/2, -dom_z/2-PML_bot, PML_lc};
Point(10)  = {-dom_x/2,  dom_y/2, -dom_z/2-PML_bot, PML_lc};
Point(11)  = { dom_x/2,  dom_y/2, -dom_z/2-PML_bot, PML_lc};
Point(12)  = { dom_x/2, -dom_y/2, -dom_z/2-PML_bot, PML_lc};
Point(13)  = {-dom_x/2, -dom_y/2,  dom_z/2+PML_top, PML_lc};
Point(14)  = {-dom_x/2,  dom_y/2,  dom_z/2+PML_top, PML_lc};
Point(15)  = { dom_x/2,  dom_y/2,  dom_z/2+PML_top, PML_lc};
Point(16)  = { dom_x/2, -dom_y/2,  dom_z/2+PML_top, PML_lc};

// PML along Y
Point(17)  = {-dom_x/2, -dom_y/2-PML_lat, -dom_z/2, PML_lc};
Point(18)  = {-dom_x/2,  dom_y/2+PML_lat, -dom_z/2, PML_lc};
Point(19)  = { dom_x/2,  dom_y/2+PML_lat, -dom_z/2, PML_lc};
Point(20)  = { dom_x/2, -dom_y/2-PML_lat, -dom_z/2, PML_lc};
Point(21)  = {-dom_x/2, -dom_y/2-PML_lat,  dom_z/2, PML_lc};
Point(22)  = {-dom_x/2,  dom_y/2+PML_lat,  dom_z/2, PML_lc};
Point(23)  = { dom_x/2,  dom_y/2+PML_lat,  dom_z/2, PML_lc};
Point(24)  = { dom_x/2, -dom_y/2-PML_lat,  dom_z/2, PML_lc};

// PML along X
Point(25)  = {-dom_x/2-PML_lat, -dom_y/2, -dom_z/2, PML_lc};
Point(26)  = {-dom_x/2-PML_lat,  dom_y/2, -dom_z/2, PML_lc};
Point(27)  = { dom_x/2+PML_lat,  dom_y/2, -dom_z/2, PML_lc};
Point(28)  = { dom_x/2+PML_lat, -dom_y/2, -dom_z/2, PML_lc};
Point(29)  = {-dom_x/2-PML_lat, -dom_y/2,  dom_z/2, PML_lc};
Point(30)  = {-dom_x/2-PML_lat,  dom_y/2,  dom_z/2, PML_lc};
Point(31)  = { dom_x/2+PML_lat,  dom_y/2,  dom_z/2, PML_lc};
Point(32)  = { dom_x/2+PML_lat, -dom_y/2,  dom_z/2, PML_lc};

// 8 missing points
Point(33)  = {-dom_x/2-PML_lat, -dom_y/2-PML_lat, -dom_z/2-PML_bot, PML_lc};
Point(34)  = {-dom_x/2-PML_lat,  dom_y/2+PML_lat, -dom_z/2-PML_bot, PML_lc};
Point(35)  = { dom_x/2+PML_lat,  dom_y/2+PML_lat, -dom_z/2-PML_bot, PML_lc};
Point(36)  = { dom_x/2+PML_lat, -dom_y/2-PML_lat, -dom_z/2-PML_bot, PML_lc};
Point(37)  = {-dom_x/2-PML_lat, -dom_y/2-PML_lat,  dom_z/2+PML_top, PML_lc};
Point(38)  = {-dom_x/2-PML_lat,  dom_y/2+PML_lat,  dom_z/2+PML_top, PML_lc};
Point(39)  = { dom_x/2+PML_lat,  dom_y/2+PML_lat,  dom_z/2+PML_top, PML_lc};
Point(40)  = { dom_x/2+PML_lat, -dom_y/2-PML_lat,  dom_z/2+PML_top, PML_lc};


Point(41)  = {-dom_x/2, -dom_y/2-PML_lat, -dom_z/2-PML_bot, PML_lc};
Point(42)  = {-dom_x/2,  dom_y/2+PML_lat, -dom_z/2-PML_bot, PML_lc};
Point(43)  = { dom_x/2, -dom_y/2-PML_lat, -dom_z/2-PML_top, PML_lc};
Point(44)  = { dom_x/2,  dom_y/2+PML_lat, -dom_z/2-PML_top, PML_lc};
Point(45)  = {-dom_x/2, -dom_y/2-PML_lat,  dom_z/2+PML_bot, PML_lc};
Point(46)  = {-dom_x/2,  dom_y/2+PML_lat,  dom_z/2+PML_bot, PML_lc};
Point(47)  = { dom_x/2, -dom_y/2-PML_lat,  dom_z/2+PML_top, PML_lc};
Point(48)  = { dom_x/2,  dom_y/2+PML_lat,  dom_z/2+PML_top, PML_lc};


Point(49)  = {-dom_x/2-PML_lat, -dom_y/2, -dom_z/2-PML_bot, PML_lc};
Point(50)  = { dom_x/2+PML_lat, -dom_y/2, -dom_z/2-PML_bot, PML_lc};
Point(51)  = {-dom_x/2-PML_lat, -dom_y/2,  dom_z/2+PML_bot, PML_lc};
Point(52)  = { dom_x/2+PML_lat, -dom_y/2,  dom_z/2+PML_bot, PML_lc};
Point(53)  = {-dom_x/2-PML_lat,  dom_y/2, -dom_z/2-PML_bot, PML_lc};
Point(54)  = { dom_x/2+PML_lat,  dom_y/2, -dom_z/2-PML_bot, PML_lc};
Point(55)  = {-dom_x/2-PML_lat,  dom_y/2,  dom_z/2+PML_bot, PML_lc};
Point(56)  = { dom_x/2+PML_lat,  dom_y/2,  dom_z/2+PML_bot, PML_lc};

Point(57)  = {-dom_x/2-PML_lat, -dom_y/2-PML_lat, -dom_z/2, PML_lc};
Point(58)  = { dom_x/2+PML_lat, -dom_y/2-PML_lat, -dom_z/2, PML_lc};
Point(59)  = {-dom_x/2-PML_lat,  dom_y/2+PML_lat, -dom_z/2, PML_lc};
Point(60)  = { dom_x/2+PML_lat,  dom_y/2+PML_lat, -dom_z/2, PML_lc};
Point(61)  = {-dom_x/2-PML_lat, -dom_y/2-PML_lat,  dom_z/2, PML_lc};
Point(62)  = { dom_x/2+PML_lat, -dom_y/2-PML_lat,  dom_z/2, PML_lc};
Point(63)  = {-dom_x/2-PML_lat,  dom_y/2+PML_lat,  dom_z/2, PML_lc};
Point(64)  = { dom_x/2+PML_lat,  dom_y/2+PML_lat,  dom_z/2, PML_lc};

// Sphere
Point(80) = { 0,  0, ro, In_lc};
Point(81) = { 0,  0,-ro, In_lc};
Point(82) = { 0, ro,  0, In_lc};
Point(83) = { 0,-ro,  0, In_lc};
Point(85) = { ro, 0,  0, In_lc};
Point(86) = {-ro, 0,  0, In_lc};
Point(87) = {  0, 0,  0, CenterScat_lc};

Line(1)   = {33, 49};
Line(2)   = {49, 53};
Line(3)   = {53, 34};
Line(4)   = {34, 42};
Line(5)   = {42, 44};
Line(6)   = {44, 35};
Line(7)   = {35, 54};
Line(8)   = {54, 11};
Line(9)   = {11, 10};
Line(10)  = {10, 53};
Line(11)  = {49,  9};
Line(12)  = { 9, 12};
Line(13)  = {12, 50};
Line(14)  = {50, 54};
Line(15)  = {50, 36};
Line(16)  = {36, 43};
Line(17)  = {43, 41};
Line(18)  = {41, 33};
Line(19)  = {41,  9};
Line(20)  = { 9, 10};
Line(21)  = {10, 42};
Line(22)  = {43, 12};
Line(23)  = {12, 11};
Line(24)  = {11, 44};
Line(25)  = {57, 25};
Line(26)  = {25, 26};
Line(27)  = {26, 59};
Line(28)  = {59, 18};
Line(29)  = {18, 19};
Line(30)  = {19, 60};
Line(31)  = {60, 27};
Line(32)  = {27, 28};
Line(33)  = {28, 58};
Line(34)  = {58, 20};
Line(35)  = {20, 17};
Line(36)  = {17, 57};
Line(37)  = {25,  1};
Line(38)  = { 1,  4};
Line(39)  = { 4, 28};
Line(40)  = {27,  3};
Line(41)  = { 3,  2};
Line(42)  = { 2, 26};
Line(43)  = {17,  1};
Line(44)  = { 1,  2};
Line(45)  = { 2, 18};
Line(46)  = {19,  3};
Line(47)  = { 3,  4};
Line(48)  = { 4, 20};
Line(49)  = {61, 21};
Line(50)  = {21, 24};
Line(51)  = {24, 62};
Line(52)  = {62, 32};
Line(53)  = {32, 31};
Line(54)  = {31, 64};
Line(55)  = {64, 23};
Line(56)  = {23, 22};
Line(57)  = {22, 63};
Line(58)  = {63, 30};
Line(59)  = {30, 29};
Line(60)  = {29, 61};
Line(61)  = {21,  5};
Line(62)  = { 5,  6};
Line(63)  = { 6, 22};
Line(64)  = {23,  7};
Line(65)  = { 7,  8};
Line(66)  = { 8, 24};
Line(67)  = {32,  8};
Line(68)  = { 8,  5};
Line(69)  = { 5, 29};
Line(70)  = {30,  6};
Line(71)  = { 6,  7};
Line(72)  = { 7, 31};
Line(73)  = {37, 45};
Line(74)  = {45, 47};
Line(75)  = {47, 40};
Line(76)  = {40, 52};
Line(77)  = {52, 56};
Line(78)  = {56, 39};
Line(79)  = {39, 48};
Line(80)  = {48, 46};
Line(81)  = {46, 38};
Line(82)  = {38, 55};
Line(83)  = {55, 51};
Line(84)  = {51, 37};
Line(85)  = {45, 13};
Line(86)  = {13, 14};
Line(87)  = {14, 46};
Line(88)  = {55, 14};
Line(89)  = {14, 15};
Line(90)  = {15, 56};
Line(91)  = {52, 16};
Line(92)  = {16, 13};
Line(93)  = {13, 51};
Line(94)  = {48, 15};
Line(95)  = {15, 16};
Line(96)  = {16, 47};
Line(97)  = {37, 61};
Line(98)  = {61, 57};
Line(99)  = {57, 33};
Line(100) = {45, 21};
Line(101) = {21, 17};
Line(102) = {17, 41};
Line(103) = {47, 24};
Line(104) = {24, 20};
Line(105) = {20, 43};
Line(106) = {40, 62};
Line(107) = {62, 58};
Line(108) = {58, 36};
Line(109) = {52, 32};
Line(110) = {32, 28};
Line(111) = {28, 50};
Line(112) = {56, 31};
Line(113) = {31, 27};
Line(114) = {27, 54};
Line(115) = {39, 64};
Line(116) = {64, 60};
Line(117) = {60, 35};

Line(118) = {48, 23};
Line(119) = {23, 19};
Line(120) = {19, 44};
Line(121) = {46, 22};
Line(122) = {22, 18};
Line(123) = {18, 42};
Line(124) = {38, 63};
Line(125) = {63, 59};
Line(126) = {59, 34};
Line(127) = {55, 30};
Line(128) = {30, 26};
Line(129) = {26, 53};
Line(130) = {51, 29};
Line(131) = {29, 25};
Line(132) = {25, 49};
Line(133) = {14,  6};
Line(134) = { 6,  2};
Line(135) = { 2, 10};
Line(136) = {15,  7};
Line(137) = { 7,  3};
Line(138) = { 3, 11};
Line(139) = {16,  8};
Line(140) = { 8,  4};
Line(141) = { 4, 12};
Line(142) = {13,  5};
Line(143) = { 5,  1};
Line(144) = { 1,  9};

Circle(145) = {82, 87, 86};
Circle(146) = {86, 87, 83};
Circle(147) = {83, 87, 85};
Circle(148) = {85, 87, 82};
Circle(149) = {81, 87, 86};
Circle(150) = {86, 87, 80};
Circle(151) = {80, 87, 85};
Circle(152) = {85, 87, 81};
Circle(153) = {83, 87, 80};
Circle(154) = {80, 87, 82};
Circle(155) = {82, 87, 81};
Circle(156) = {81, 87, 83};
Line Loop(157) = {151, -147, 153};
Ruled Surface(158) = {157};
Line Loop(159) = {151, 148, -154};
Ruled Surface(160) = {159};
Line Loop(161) = {145, 150, 154};
Ruled Surface(162) = {161};
Line Loop(163) = {150, -153, -146};
Ruled Surface(164) = {163};
Line Loop(165) = {156, 147, 152};
Ruled Surface(166) = {165};
Line Loop(167) = {152, -155, -148};
Ruled Surface(168) = {167};
Line Loop(169) = {155, 149, -145};
Ruled Surface(170) = {169};
Line Loop(171) = {149, 146, -156};
Ruled Surface(172) = {171};
Line Loop(173) = {80, -87, 89, -94};
Plane Surface(174) = {173};
Line Loop(175) = {94, 90, 78, 79};
Plane Surface(176) = {175};
Line Loop(177) = {77, -90, 95, -91};
Plane Surface(178) = {177};
Line Loop(179) = {95, 92, 86, 89};
Plane Surface(180) = {179};
Line Loop(181) = {91, 96, 75, 76};
Plane Surface(182) = {181};
Line Loop(183) = {92, -85, 74, -96};
Plane Surface(184) = {183};
Line Loop(185) = {82, 88, 87, 81};
Plane Surface(186) = {185};
Line Loop(187) = {86, -88, 83, -93};
Plane Surface(188) = {187};
Line Loop(189) = {93, 84, 73, 85};
Plane Surface(190) = {189};
Line Loop(191) = {58, 70, 63, 57};
Plane Surface(192) = {191};
Line Loop(193) = {63, -56, 64, -71};
Plane Surface(194) = {193};
Line Loop(195) = {64, 72, 54, 55};
Plane Surface(196) = {195};
Line Loop(197) = {72, -53, 67, -65};
Plane Surface(198) = {197};
Line Loop(199) = {67, 66, 51, 52};
Plane Surface(200) = {199};
Line Loop(201) = {71, 65, 68, 62};
Plane Surface(202) = {201};
Line Loop(203) = {68, -61, 50, -66};
Plane Surface(204) = {203};
Line Loop(205) = {59, -69, 62, -70};
Plane Surface(206) = {205};
Line Loop(207) = {69, 60, 49, 61};
Plane Surface(208) = {207};
Line Loop(209) = {28, -45, 42, 27};
Plane Surface(210) = {209};
Line Loop(211) = {29, 46, 41, 45};
Plane Surface(212) = {211};
Line Loop(213) = {30, 31, 40, -46};
Plane Surface(214) = {213};
Line Loop(215) = {40, 47, 39, -32};
Plane Surface(216) = {215};
Line Loop(217) = {41, -44, 38, -47};
Plane Surface(218) = {217};
Line Loop(219) = {39, 33, 34, -48};
Plane Surface(220) = {219};
Line Loop(221) = {38, 48, 35, 43};
Plane Surface(222) = {221};
Line Loop(223) = {26, -42, -44, -37};
Plane Surface(224) = {223};
Line Loop(225) = {37, -43, 36, 25};
Plane Surface(226) = {225};
Line Loop(227) = {3, 4, -21, 10};
Plane Surface(228) = {227};
Line Loop(229) = {5, -24, 9, 21};
Plane Surface(230) = {229};
Line Loop(231) = {6, 7, 8, 24};
Plane Surface(232) = {231};
Line Loop(233) = {14, 8, -23, 13};
Plane Surface(234) = {233};
Line Loop(235) = {9, -20, 12, 23};
Plane Surface(236) = {235};
Line Loop(237) = {2, -10, -20, -11};
Plane Surface(238) = {237};
Line Loop(239) = {13, 15, 16, 22};
Plane Surface(240) = {239};
Line Loop(241) = {12, -22, 17, 19};
Plane Surface(242) = {241};
Line Loop(243) = {19, -11, -1, -18};
Plane Surface(244) = {243};
Line Loop(245) = {127, 70, -133, -88};
Plane Surface(246) = {245};
Line Loop(247) = {133, 71, -136, -89};
Plane Surface(248) = {247};
Line Loop(249) = {90, 112, -72, -136};
Plane Surface(250) = {249};
Line Loop(251) = {128, -42, -134, -70};
Plane Surface(252) = {251};
Line Loop(253) = {41, -134, 71, 137};
Plane Surface(254) = {253};
Line Loop(255) = {137, -40, -113, -72};
Plane Surface(256) = {255};
Line Loop(257) = {129, -10, -135, 42};
Plane Surface(258) = {257};
Line Loop(259) = {135, -9, -138, 41};
Plane Surface(260) = {259};
Line Loop(261) = {138, -8, -114, 40};
Plane Surface(262) = {261};
Line Loop(263) = {130, -69, -142, 93};
Plane Surface(264) = {263};
Line Loop(265) = {92, 142, -68, -139};
Plane Surface(266) = {265};
Line Loop(267) = {139, -67, -109, 91};
Plane Surface(268) = {267};
Line Loop(269) = {131, 37, -143, 69};
Plane Surface(270) = {269};
Line Loop(271) = {143, 38, -140, 68};
Plane Surface(272) = {271};
Line Loop(273) = {140, 39, -110, 67};
Plane Surface(274) = {273};
Line Loop(275) = {132, 11, -144, -37};
Plane Surface(276) = {275};
Line Loop(277) = {12, -141, -38, 144};
Plane Surface(278) = {277};
Line Loop(279) = {141, 13, -111, -39};
Plane Surface(280) = {279};
Line Loop(281) = {99, -18, -102, 36};
Plane Surface(282) = {281};
Line Loop(283) = {102, -17, -105, 35};
Plane Surface(284) = {283};
Line Loop(285) = {34, 105, -16, -108};
Plane Surface(286) = {285};
Line Loop(287) = {34, -104, 51, 107};
Plane Surface(288) = {287};
Line Loop(289) = {50, 104, 35, -101};
Plane Surface(290) = {289};
Line Loop(291) = {36, -98, 49, 101};
Plane Surface(292) = {291};
Line Loop(293) = {97, 49, -100, -73};
Plane Surface(294) = {293};
Line Loop(295) = {74, 103, -50, -100};
Plane Surface(296) = {295};
Line Loop(297) = {103, 51, -106, -75};
Plane Surface(298) = {297};
Line Loop(299) = {126, 4, -123, -28};
Plane Surface(300) = {299};
Line Loop(301) = {125, 28, -122, 57};
Plane Surface(302) = {301};
Line Loop(303) = {122, 29, -119, 56};
Plane Surface(304) = {303};
Line Loop(305) = {119, 30, -116, 55};
Plane Surface(306) = {305};
Line Loop(307) = {55, -118, -79, 115};
Plane Surface(308) = {307};
Line Loop(309) = {80, 121, -56, -118};
Plane Surface(310) = {309};
Line Loop(311) = {81, 124, -57, -121};
Plane Surface(312) = {311};
Line Loop(313) = {123, 5, -120, -29};
Plane Surface(314) = {313};
Line Loop(315) = {120, 6, -117, -30};
Plane Surface(316) = {315};
Line Loop(317) = {121, -63, -133, 87};
Plane Surface(318) = {317};
Line Loop(319) = {133, -62, -142, 86};
Plane Surface(320) = {319};
Line Loop(321) = {142, -61, -100, 85};
Plane Surface(322) = {321};
Line Loop(323) = {122, -45, -134, 63};
Plane Surface(324) = {323};
Line Loop(325) = {134, -44, -143, 62};
Plane Surface(326) = {325};
Line Loop(327) = {143, -43, -101, 61};
Plane Surface(328) = {327};
Line Loop(329) = {123, -21, -135, 45};
Plane Surface(330) = {329};
Line Loop(331) = {135, -20, -144, 44};
Plane Surface(332) = {331};
Line Loop(333) = {144, -19, -102, 43};
Plane Surface(334) = {333};
Line Loop(335) = {126, -3, -129, 27};
Plane Surface(336) = {335};
Line Loop(337) = {129, -2, -132, 26};
Plane Surface(338) = {337};
Line Loop(339) = {132, -1, -99, 25};
Plane Surface(340) = {339};
Line Loop(341) = {125, -27, -128, -58};
Plane Surface(342) = {341};
Line Loop(343) = {128, -26, -131, -59};
Plane Surface(344) = {343};
Line Loop(345) = {131, -25, -98, -60};
Plane Surface(346) = {345};
Line Loop(347) = {97, -60, -130, 84};
Plane Surface(348) = {347};
Line Loop(349) = {83, 130, -59, -127};
Plane Surface(350) = {349};
Line Loop(351) = {58, -127, -82, 124};
Plane Surface(352) = {351};
Line Loop(353) = {96, 103, -66, -139};
Plane Surface(354) = {353};
Line Loop(355) = {139, -65, -136, 95};
Plane Surface(356) = {355};
Line Loop(357) = {136, -64, -118, 94};
Plane Surface(358) = {357};
Line Loop(359) = {66, 104, -48, -140};
Plane Surface(360) = {359};
Line Loop(361) = {140, -47, -137, 65};
Plane Surface(362) = {361};
Line Loop(363) = {137, -46, -119, 64};
Plane Surface(364) = {363};
Line Loop(365) = {105, 22, -141, 48};
Plane Surface(366) = {365};
Line Loop(367) = {141, 23, -138, 47};
Plane Surface(368) = {367};
Line Loop(369) = {138, 24, -120, 46};
Plane Surface(370) = {369};
Line Loop(371) = {108, -15, -111, 33};
Plane Surface(372) = {371};
Line Loop(373) = {111, 14, -114, 32};
Plane Surface(374) = {373};
Line Loop(375) = {114, -7, -117, 31};
Plane Surface(376) = {375};
Line Loop(377) = {107, -33, -110, -52};
Plane Surface(378) = {377};
Line Loop(379) = {110, -32, -113, -53};
Plane Surface(380) = {379};
Line Loop(381) = {113, -31, -116, -54};
Plane Surface(382) = {381};
Line Loop(383) = {106, 52, -109, -76};
Plane Surface(384) = {383};
Line Loop(385) = {109, 53, -112, -77};
Plane Surface(386) = {385};
Line Loop(387) = {112, 54, -115, -78};
Plane Surface(388) = {387};

Surface Loop(389) = {164, 162, 170, 168, 166, 172, 158, 160};
Volume(390) = {389};
Surface Loop(391) = {218, 272, 362, 254, 326, 202};
Volume(392) = {389, 391};
Surface Loop(393) = {182, 298, 384, 354, 200, 268};
Volume(394) = {393};
Surface Loop(395) = {184, 296, 354, 266, 204, 322};
Volume(396) = {395};
Surface Loop(397) = {190, 348, 294, 264, 208, 322};
Volume(398) = {397};
Surface Loop(399) = {188, 350, 264, 246, 320, 206};
Volume(400) = {399};
Surface Loop(401) = {186, 352, 312, 192, 318, 246};
Volume(402) = {401};
Surface Loop(403) = {180, 248, 266, 202, 356, 320};
Volume(404) = {403};
Surface Loop(405) = {310, 174, 318, 248, 194, 358};
Volume(406) = {405};
Surface Loop(407) = {176, 388, 308, 358, 250, 196};
Volume(408) = {407};
Surface Loop(409) = {178, 386, 198, 356, 268, 250};
Volume(410) = {409};
Surface Loop(411) = {378, 288, 200, 274, 360, 220};
Volume(412) = {411};
Surface Loop(413) = {216, 380, 198, 274, 256, 362};
Volume(414) = {413};
Surface Loop(415) = {306, 382, 214, 364, 196, 256};
Volume(416) = {415};
Surface Loop(417) = {304, 194, 364, 212, 254, 324};
Volume(418) = {417};
Surface Loop(419) = {302, 342, 192, 324, 252, 210};
Volume(420) = {419};
Surface Loop(421) = {344, 206, 252, 224, 270, 326};
Volume(422) = {421};
Surface Loop(423) = {292, 346, 328, 208, 270, 226};
Volume(424) = {423};
Surface Loop(425) = {290, 360, 204, 272, 328, 222};
Volume(426) = {425};
Surface Loop(427) = {240, 372, 286, 220, 366, 280};
Volume(428) = {427};
Surface Loop(429) = {234, 374, 216, 368, 262, 280};
Volume(430) = {429};
Surface Loop(431) = {376, 232, 316, 370, 214, 262};
Volume(432) = {431};
Surface Loop(433) = {314, 230, 260, 212, 330, 370};
Volume(434) = {433};
Surface Loop(435) = {228, 336, 300, 210, 258, 330};
Volume(436) = {435};
Surface Loop(437) = {236, 332, 368, 278, 260, 218};
Volume(438) = {437};
Surface Loop(439) = {284, 242, 366, 334, 222, 278};
Volume(440) = {439};
Surface Loop(441) = {282, 340, 244, 334, 226, 276};
Volume(442) = {441};
Surface Loop(443) = {238, 338, 258, 276, 332, 224};
Volume(444) = {443};

Physical Volume(1000) = {402, 398, 394, 408, 442, 428, 432, 436}; // PML XYZ
Physical Volume(1001) = {400, 410, 430, 444}; // PML XZ
Physical Volume(1002) = {406, 396, 434, 440}; // PML YZ
Physical Volume(1003) = {412, 416, 420, 424}; // PML XY
Physical Volume(1004) = {404, 438}; // PML Z
Physical Volume(1005) = {426, 418}; // PML Y
Physical Volume(1006) = {414, 422}; // PML X
Physical Volume(1007) = {392}; // Out
Physical Volume(1008) = {390}; // Scatterer

Mesh.Algorithm   = 1; // // 1=MeshAdapt, 5=Delaunay, 6=Frontal
Mesh.Algorithm3D = 1;//4; // // 1=Delaunay, 4=Frontal
Mesh.Optimize = 1;
