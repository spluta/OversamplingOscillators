#include <array>

std::array<float, 1000> getSergeWavetable() {
    std::array<float, 1000> arr = {-1.0, -0.97936677932739, -0.95382595062256, -0.92388689517975, -0.89003127813339, -0.85271441936493, -0.81236624717712, -0.76939231157303, -0.72417414188385, -0.67707097530365, -0.62842011451721, -0.578537940979, -0.52772086858749, -0.47624570131302, -0.42437115311623, -0.37233778834343, -0.32036954164505, -0.26867371797562, -0.21744230389595, -0.16685220599174, -0.1170661225915, -0.068233236670494, -0.020489417016506, 0.02604166045785, 0.071248166263103, 0.11502959579229, 0.15729621052742, 0.197968557477, 0.23697650432587, 0.27425986528397, 0.3097665309906, 0.34345278143883, 0.3752832710743, 0.40522947907448, 0.43327024579048, 0.45939087867737, 0.4835829436779, 0.50584369897842, 0.52617639303207, 0.544588804245, 0.56109398603439, 0.5757092833519, 0.58845633268356, 0.599360704422, 0.60845142602921, 0.61576104164124, 0.62132501602173, 0.62518185377121, 0.62737244367599, 0.62794011831284, 0.62693029642105, 0.62439042329788, 0.62036955356598, 0.6149183511734, 0.6080886721611, 0.59993356466293, 0.59050738811493, 0.5798647403717, 0.56806129217148, 0.55515319108963, 0.54119700193405, 0.52624934911728, 0.51036703586578, 0.493606954813, 0.47602587938309, 0.45768025517464, 0.43862646818161, 0.4189201593399, 0.39861670136452, 0.37777072191238, 0.35643631219864, 0.33466675877571, 0.31251457333565, 0.29003125429153, 0.26726767420769, 0.24427342414856, 0.22109703719616, 0.19778627157211, 0.17438752949238, 0.15094594657421, 0.12750570476055, 0.10410954803228, 0.080799154937267, 0.057614635676146, 0.034595031291246, 0.011777936480939, -0.010800368152559, -0.033104937523603, -0.055102337151766, -0.076760306954384, -0.098048143088818, -0.11893640458584, -0.13939723372459, -0.15940390527248, -0.17893129587173, -0.19795551896095, -0.21645423769951, -0.23440629243851, -0.25179201364517, -0.26859298348427, -0.28479227423668, -0.30037409067154, -0.31532406806946, -0.32962906360626, -0.34327730536461, -0.35625803470612, -0.36856198310852, -0.38018101453781, -0.39110821485519, -0.40133759379387, -0.41086459159851, -0.41968566179276, -0.42779830098152, -0.43520113825798, -0.44189381599426, -0.44787690043449, -0.4531521499157, -0.4577220082283, -0.46159017086029, -0.46476086974144, -0.4672394990921, -0.46903213858604, -0.47014567255974, -0.47058805823326, -0.47036749124527, -0.46949335932732, -0.46797552704811, -0.46582448482513, -0.46305149793625, -0.45966848731041, -0.45568764209747, -0.45112210512161, -0.44598519802094, -0.44029101729393, -0.43405395746231, -0.4272888302803, -0.42001089453697, -0.41223594546318, -0.40397995710373, -0.39525917172432, -0.38609033823013, -0.37649032473564, -0.36647635698318, -0.35606580972672, -0.34527623653412, -0.33412551879883, -0.32263135910034, -0.31081208586693, -0.29868566989899, -0.28627038002014, -0.27358451485634, -0.2606463432312, -0.24747431278229, -0.23408664762974, -0.22050173580647, -0.20673783123493, -0.19281309843063, -0.17874559760094, -0.16455344855785, -0.15025442838669, -0.13586632907391, -0.12140678614378, -0.1068931594491, -0.092342734336853, -0.077772423624992, -0.06319922208786, -0.04863964766264, -0.034109935164452, -0.019626345485449, -0.0052045620977879, 0.0091398051008582, 0.023391615599394, 0.037535730749369, 0.051557775586843, 0.065443247556686, 0.079178243875504, 0.092749170958996, 0.10614275932312, 0.11934601515532, 0.13234642148018, 0.14513179659843, 0.157690346241, 0.17001059651375, 0.18208165466785, 0.19389268755913, 0.20543356239796, 0.21669441461563, 0.22766584157944, 0.23833879828453, 0.24870455265045, 0.25875517725945, 0.26848259568214, 0.27787965536118, 0.28693923354149, 0.29565498232841, 0.30402058362961, 0.31203046441078, 0.31967929005623, 0.32696217298508, 0.33387476205826, 0.34041291475296, 0.34657302498817, 0.35235190391541, 0.35774672031403, 0.3627550303936, 0.36737489700317, 0.37160462141037, 0.37544298171997, 0.37888917326927, 0.38194274902344, 0.38460364937782, 0.38687211275101, 0.38874879479408, 0.39023479819298, 0.39133143424988, 0.39204052090645, 0.39236405491829, 0.39230453968048, 0.39186468720436, 0.39104756712914, 0.38985660672188, 0.3882954120636, 0.38636818528175, 0.38407909870148, 0.38143280148506, 0.37843412160873, 0.37508833408356, 0.37140080332756, 0.36737713217735, 0.36302343010902, 0.35834577679634, 0.35335057973862, 0.34804448485374, 0.34243443608284, 0.33652740716934, 0.33033064007759, 0.32385167479515, 0.31709814071655, 0.31007793545723, 0.30279892683029, 0.2952693104744, 0.28749740123749, 0.27949160337448, 0.27126055955887, 0.26281288266182, 0.25415745377541, 0.24530313909054, 0.23625910282135, 0.22703416645527, 0.21763777732849, 0.20807914435863, 0.19836749136448, 0.18851232528687, 0.17852300405502, 0.16840904951096, 0.15817986428738, 0.14784508943558, 0.13741421699524, 0.12689681351185, 0.11630239337683, 0.10564056038857, 0.094920761883259, 0.084152571856976, 0.073345422744751, 0.062508769333363, 0.05165196210146, 0.040784355252981, 0.0299152713269, 0.019053952768445, 0.0082094054669142, -0.0026091914623976, -0.013392901979387, -0.024132834747434, -0.034820217639208, -0.04544635489583, -0.056002613157034, -0.066480688750744, -0.076872043311596, -0.087168626487255, -0.097362294793129, -0.10744506120682, -0.11740925163031, -0.12724709510803, -0.13695111870766, -0.14651396870613, -0.15592853724957, -0.16518767178059, -0.17428462207317, -0.18321262300014, -0.19196528196335, -0.20053619146347, -0.20891919732094, -0.21710839867592, -0.22509796917439, -0.2328823953867, -0.24045622348785, -0.24781431257725, -0.25495171546936, -0.26186355948448, -0.26854529976845, -0.27499261498451, -0.2812012732029, -0.28716737031937, -0.29288718104362, -0.29835703969002, -0.30357372760773, -0.30853414535522, -0.31323537230492, -0.31767463684082, -0.32184967398643, -0.32575815916061, -0.32939800620079, -0.33276736736298, -0.33586478233337, -0.33868870139122, -0.34123814105988, -0.34351196885109, -0.34550958871841, -0.34723034501076, -0.34867402911186, -0.34984046220779, -0.35072979331017, -0.3513423204422, -0.35167846083641, -0.35173913836479, -0.35152506828308, -0.35103759169579, -0.35027784109116, -0.34924739599228, -0.34794801473618, -0.34638157486916, -0.34455010294914, -0.34245604276657, -0.34010165929794, -0.33748978376389, -0.33462309837341, -0.33150473237038, -0.32813769578934, -0.32452550530434, -0.32067158818245, -0.31657958030701, -0.31225326657295, -0.30769675970078, -0.302914083004, -0.29790955781937, -0.29268765449524, -0.28725284337997, -0.28160980343819, -0.27576351165771, -0.26971879601479, -0.26348075270653, -0.25705462694168, -0.25044569373131, -0.24365937709808, -0.23670123517513, -0.22957690060139, -0.22229206562042, -0.21485260128975, -0.2072644084692, -0.19953340291977, -0.191665828228, -0.18366770446301, -0.17554526031017, -0.1673047542572, -0.1589527130127, -0.15049533545971, -0.14193916320801, -0.13329072296619, -0.12455654889345, -0.11574317514896, -0.10685729235411, -0.097905471920967, -0.088894419372082, -0.079830840229988, -0.070721432566643, -0.061572916805744, -0.052391976118088, -0.04318530485034, -0.033959724009037, -0.02472186461091, -0.015478426590562, -0.0062361257150769, 0.0029984398279339, 0.012218607589602, 0.021417772397399, 0.030589360743761, 0.039726927876472, 0.048823848366737, 0.057873856276274, 0.066870398819447, 0.075807340443134, 0.084678180515766, 0.093476921319962, 0.10219732671976, 0.11083338409662, 0.11937906593084, 0.12782843410969, 0.13617573678493, 0.14441521465778, 0.15254113078117, 0.16054801642895, 0.16843041777611, 0.17618297040462, 0.18380033969879, 0.19127748906612, 0.19860927760601, 0.20579083263874, 0.21281731128693, 0.21968403458595, 0.22638639807701, 0.23292005062103, 0.23928053677082, 0.24546374380589, 0.25146555900574, 0.25728207826614, 0.262909501791, 0.26834419369698, 0.27358260750771, 0.27862137556076, 0.28345721960068, 0.28808709979057, 0.29250809550285, 0.29671737551689, 0.30071228742599, 0.30449035763741, 0.30804923176765, 0.31138679385185, 0.31450092792511, 0.31738975644112, 0.32005164027214, 0.32248491048813, 0.32468819618225, 0.32666033506393, 0.32840013504028, 0.32990676164627, 0.33117932081223, 0.33221733570099, 0.33302024006844, 0.33358785510063, 0.33392003178596, 0.33401668071747, 0.33387809991837, 0.3335046172142, 0.33289670944214, 0.33205506205559, 0.33098050951958, 0.32967403531075, 0.3281367123127, 0.32636985182762, 0.32437488436699, 0.32215335965157, 0.31970712542534, 0.31703796982765, 0.31414800882339, 0.31103929877281, 0.30771428346634, 0.30417537689209, 0.30042517185211, 0.29646646976471, 0.29230210185051, 0.28793516755104, 0.28336870670319, 0.27860605716705, 0.27365064620972, 0.2685059607029, 0.26317572593689, 0.257663667202, 0.25197371840477, 0.24610990285873, 0.24007628858089, 0.23387718200684, 0.22751688957214, 0.22099989652634, 0.21433073282242, 0.20751401782036, 0.20055454969406, 0.19345712661743, 0.18622675538063, 0.17886836826801, 0.17138704657555, 0.16378806531429, 0.15607658028603, 0.14825797080994, 0.14033761620522, 0.13232097029686, 0.12421354651451, 0.11602096259594, 0.10774881392717, 0.09940280020237, 0.090988650918007, 0.082512147724628, 0.073979109525681, 0.065395385026932, 0.056766856461763, 0.048099439591169, 0.039399087429047, 0.030671752989292, 0.021923407912254, 0.013160053640604, 0.0043876855634153, -0.0043876855634153, -0.013160053640604, -0.021923407912254, -0.030671752989292, -0.039399087429047, -0.048099439591169, -0.056766856461763, -0.065395385026932, -0.073979109525681, -0.082512147724628, -0.090988650918007, -0.09940280020237, -0.10774881392717, -0.11602096259594, -0.12421354651451, -0.13232097029686, -0.14033761620522, -0.14825797080994, -0.15607658028603, -0.16378806531429, -0.17138704657555, -0.17886836826801, -0.18622675538063, -0.19345712661743, -0.20055454969406, -0.20751401782036, -0.21433073282242, -0.22099989652634, -0.22751688957214, -0.23387718200684, -0.24007628858089, -0.24610990285873, -0.25197371840477, -0.257663667202, -0.26317572593689, -0.2685059607029, -0.27365064620972, -0.27860605716705, -0.28336870670319, -0.28793516755104, -0.29230210185051, -0.29646646976471, -0.30042517185211, -0.30417537689209, -0.30771428346634, -0.31103929877281, -0.31414800882339, -0.31703796982765, -0.31970712542534, -0.32215335965157, -0.32437488436699, -0.32636985182762, -0.3281367123127, -0.32967403531075, -0.33098050951958, -0.33205506205559, -0.33289670944214, -0.3335046172142, -0.33387809991837, -0.33401668071747, -0.33392003178596, -0.33358785510063, -0.33302024006844, -0.33221733570099, -0.33117932081223, -0.32990676164627, -0.32840013504028, -0.32666033506393, -0.32468819618225, -0.32248491048813, -0.32005164027214, -0.31738975644112, -0.31450092792511, -0.31138679385185, -0.30804923176765, -0.30449035763741, -0.30071228742599, -0.29671737551689, -0.29250809550285, -0.28808709979057, -0.28345721960068, -0.27862137556076, -0.27358260750771, -0.26834419369698, -0.262909501791, -0.25728207826614, -0.25146555900574, -0.24546374380589, -0.23928053677082, -0.23292005062103, -0.22638639807701, -0.21968403458595, -0.21281731128693, -0.20579083263874, -0.19860927760601, -0.19127748906612, -0.18380033969879, -0.17618297040462, -0.16843041777611, -0.16054801642895, -0.15254113078117, -0.14441521465778, -0.13617573678493, -0.12782843410969, -0.11937906593084, -0.11083338409662, -0.10219732671976, -0.093476921319962, -0.084678180515766, -0.075807340443134, -0.066870398819447, -0.057873856276274, -0.048823848366737, -0.039726927876472, -0.030589360743761, -0.021417772397399, -0.012218607589602, -0.0029984398279339, 0.0062361257150769, 0.015478426590562, 0.02472186461091, 0.033959724009037, 0.04318530485034, 0.052391976118088, 0.061572916805744, 0.070721432566643, 0.079830840229988, 0.088894419372082, 0.097905471920967, 0.10685729235411, 0.11574317514896, 0.12455654889345, 0.13329072296619, 0.14193916320801, 0.15049533545971, 0.1589527130127, 0.1673047542572, 0.17554526031017, 0.18366770446301, 0.191665828228, 0.19953340291977, 0.2072644084692, 0.21485260128975, 0.22229206562042, 0.22957690060139, 0.23670123517513, 0.24365937709808, 0.25044569373131, 0.25705462694168, 0.26348075270653, 0.26971879601479, 0.27576351165771, 0.28160980343819, 0.28725284337997, 0.29268765449524, 0.29790955781937, 0.302914083004, 0.30769675970078, 0.31225326657295, 0.31657958030701, 0.32067158818245, 0.32452550530434, 0.32813769578934, 0.33150473237038, 0.33462309837341, 0.33748978376389, 0.34010165929794, 0.34245604276657, 0.34455010294914, 0.34638157486916, 0.34794801473618, 0.34924739599228, 0.35027784109116, 0.35103759169579, 0.35152506828308, 0.35173913836479, 0.35167846083641, 0.3513423204422, 0.35072979331017, 0.34984046220779, 0.34867402911186, 0.34723034501076, 0.34550958871841, 0.34351196885109, 0.34123814105988, 0.33868870139122, 0.33586478233337, 0.33276736736298, 0.32939800620079, 0.32575815916061, 0.32184967398643, 0.31767463684082, 0.31323537230492, 0.30853414535522, 0.30357372760773, 0.29835703969002, 0.29288718104362, 0.28716737031937, 0.2812012732029, 0.27499261498451, 0.26854529976845, 0.26186355948448, 0.25495171546936, 0.24781431257725, 0.24045622348785, 0.2328823953867, 0.22509796917439, 0.21710839867592, 0.20891919732094, 0.20053619146347, 0.19196528196335, 0.18321262300014, 0.17428462207317, 0.16518767178059, 0.15592853724957, 0.14651396870613, 0.13695111870766, 0.12724709510803, 0.11740925163031, 0.10744506120682, 0.097362294793129, 0.087168626487255, 0.076872043311596, 0.066480688750744, 0.056002613157034, 0.04544635489583, 0.034820217639208, 0.024132834747434, 0.013392901979387, 0.0026091914623976, -0.0082094054669142, -0.019053952768445, -0.0299152713269, -0.040784355252981, -0.05165196210146, -0.062508769333363, -0.073345422744751, -0.084152571856976, -0.094920761883259, -0.10564056038857, -0.11630239337683, -0.12689681351185, -0.13741421699524, -0.14784508943558, -0.15817986428738, -0.16840904951096, -0.17852300405502, -0.18851232528687, -0.19836749136448, -0.20807914435863, -0.21763777732849, -0.22703416645527, -0.23625910282135, -0.24530313909054, -0.25415745377541, -0.26281288266182, -0.27126055955887, -0.27949160337448, -0.28749740123749, -0.2952693104744, -0.30279892683029, -0.31007793545723, -0.31709814071655, -0.32385167479515, -0.33033064007759, -0.33652740716934, -0.34243443608284, -0.34804448485374, -0.35335057973862, -0.35834577679634, -0.36302343010902, -0.36737713217735, -0.37140080332756, -0.37508833408356, -0.37843412160873, -0.38143280148506, -0.38407909870148, -0.38636818528175, -0.3882954120636, -0.38985660672188, -0.39104756712914, -0.39186468720436, -0.39230453968048, -0.39236405491829, -0.39204052090645, -0.39133143424988, -0.39023479819298, -0.38874879479408, -0.38687211275101, -0.38460364937782, -0.38194274902344, -0.37888917326927, -0.37544298171997, -0.37160462141037, -0.36737489700317, -0.3627550303936, -0.35774672031403, -0.35235190391541, -0.34657302498817, -0.34041291475296, -0.33387476205826, -0.32696217298508, -0.31967929005623, -0.31203046441078, -0.30402058362961, -0.29565498232841, -0.28693923354149, -0.27787965536118, -0.26848259568214, -0.25875517725945, -0.24870455265045, -0.23833879828453, -0.22766584157944, -0.21669441461563, -0.20543356239796, -0.19389268755913, -0.18208165466785, -0.17001059651375, -0.157690346241, -0.14513179659843, -0.13234642148018, -0.11934601515532, -0.10614275932312, -0.092749170958996, -0.079178243875504, -0.065443247556686, -0.051557775586843, -0.037535730749369, -0.023391615599394, -0.0091398051008582, 0.0052045620977879, 0.019626345485449, 0.034109935164452, 0.04863964766264, 0.06319922208786, 0.077772423624992, 0.092342734336853, 0.1068931594491, 0.12140678614378, 0.13586632907391, 0.15025442838669, 0.16455344855785, 0.17874559760094, 0.19281309843063, 0.20673783123493, 0.22050173580647, 0.23408664762974, 0.24747431278229, 0.2606463432312, 0.27358451485634, 0.28627038002014, 0.29868566989899, 0.31081208586693, 0.32263135910034, 0.33412551879883, 0.34527623653412, 0.35606580972672, 0.36647635698318, 0.37649032473564, 0.38609033823013, 0.39525917172432, 0.40397995710373, 0.41223594546318, 0.42001089453697, 0.4272888302803, 0.43405395746231, 0.44029101729393, 0.44598519802094, 0.45112210512161, 0.45568764209747, 0.45966848731041, 0.46305149793625, 0.46582448482513, 0.46797552704811, 0.46949335932732, 0.47036749124527, 0.47058805823326, 0.47014567255974, 0.46903213858604, 0.4672394990921, 0.46476086974144, 0.46159017086029, 0.4577220082283, 0.4531521499157, 0.44787690043449, 0.44189381599426, 0.43520113825798, 0.42779830098152, 0.41968566179276, 0.41086459159851, 0.40133759379387, 0.39110821485519, 0.38018101453781, 0.36856198310852, 0.35625803470612, 0.34327730536461, 0.32962906360626, 0.31532406806946, 0.30037409067154, 0.28479227423668, 0.26859298348427, 0.25179201364517, 0.23440629243851, 0.21645423769951, 0.19795551896095, 0.17893129587173, 0.15940390527248, 0.13939723372459, 0.11893640458584, 0.098048143088818, 0.076760306954384, 0.055102337151766, 0.033104937523603, 0.010800368152559, -0.011777936480939, -0.034595031291246, -0.057614635676146, -0.080799154937267, -0.10410954803228, -0.12750570476055, -0.15094594657421, -0.17438752949238, -0.19778627157211, -0.22109703719616, -0.24427342414856, -0.26726767420769, -0.29003125429153, -0.31251457333565, -0.33466675877571, -0.35643631219864, -0.37777072191238, -0.39861670136452, -0.4189201593399, -0.43862646818161, -0.45768025517464, -0.47602587938309, -0.493606954813, -0.51036703586578, -0.52624934911728, -0.54119700193405, -0.55515319108963, -0.56806129217148, -0.5798647403717, -0.59050738811493, -0.59993356466293, -0.6080886721611, -0.6149183511734, -0.62036955356598, -0.62439042329788, -0.62693029642105, -0.62794011831284, -0.62737244367599, -0.62518185377121, -0.62132501602173, -0.61576104164124, -0.60845142602921, -0.599360704422, -0.58845633268356, -0.5757092833519, -0.56109398603439, -0.544588804245, -0.52617639303207, -0.50584369897842, -0.4835829436779, -0.45939087867737, -0.43327024579048, -0.40522947907448, -0.3752832710743, -0.34345278143883, -0.3097665309906, -0.27425986528397, -0.23697650432587, -0.197968557477, -0.15729621052742, -0.11502959579229, -0.071248166263103, -0.02604166045785, 0.020489417016506, 0.068233236670494, 0.1170661225915, 0.16685220599174, 0.21744230389595, 0.26867371797562, 0.32036954164505, 0.37233778834343, 0.42437115311623, 0.47624570131302, 0.52772086858749, 0.578537940979, 0.62842011451721, 0.67707097530365, 0.72417414188385, 0.76939231157303, 0.81236624717712, 0.85271441936493, 0.89003127813339, 0.92388689517975, 0.95382595062256, 0.97936677932739, 1.0};
    return arr;
}