# PlanetaryEphemeris abbreviation
const PE = PlanetaryEphemeris

# Integration parameters
const order = 30
const abstol = 1.0E-20

# Important bodies indexes 
const su = 1 # Sun's index
const ea = 4 # Earth's index
const mo = 5 # Moon's index

# Mass parameters μ = G*m (in au^3/day^2) of major bodies (Sun + Planets + Moon)
# See Table 8 in page 49 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const GM1 = 0.491248045036476000E-10 # Mercury
const GM2 = 0.724345233264412000E-09 # Venus
const GM3 = 0.888769244512563400E-09 # Earth
const GM4 = 0.954954869555077000E-10 # Mars
const GM5 = 0.282534584083387000E-06 # Jupiter
const GM6 = 0.845970607324503000E-07 # Saturn
const GM7 = 0.129202482578296000E-07 # Uranus
const GM8 = 0.152435734788511000E-07 # Neptune
const GM9 = 0.217844105197418000E-11 # Pluto
const GMS = 0.295912208285591100E-03 # Sun
const GMM = 0.109318945074237400E-10 # Moon
const GMB = 0.899701139019987100E-09 # Earth-Moon barycenter

# Vector of mass parameters μ = G*m (in au^3/day^2) of  Sun + Planets + Moon + 343 Asteroids
# See Table 8 in page 49 and Table 12 in pages 53-59 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const μ = [GMS, GM1, GM2, GM3, GMM, GM4, GM5, GM6, GM7, GM8, GM9,

1.4004765561723440E-13, # 1 Ceres
3.8547501878088100E-14, # 4 Vesta
3.1044481989387130E-14, # 2 Pallas
1.2358007872941250E-14, # 10 Hygiea
6.3432804736486020E-15, # 31 Euphrosyne
5.2561686784936620E-15, # 704 Interamnia
5.1981269794574980E-15, # 511 Davida
4.6783074183509050E-15, # 15 Eunomia
3.6175383171479370E-15, # 3 Juno
3.4115868261938120E-15, # 16 Psyche
3.1806592826525410E-15, # 65 Cybele
2.5771141273110470E-15, # 88 Thisbe
2.5310917260150680E-15, # 48 Doris
2.4767881012558670E-15, # 52 Europa
2.2955593906374620E-15, # 451 Patientia
2.1992951735740730E-15, # 87 Sylvia
2.1364344425714070E-15, # 7 Iris
2.1124383605999520E-15, # 423 Diotima
1.9758423651245200E-15, # 29 Amphitrite
1.8939016675253820E-15, # 24 Themis
1.7970048945074460E-15, # 13 Egeria
1.7558991833247000E-15, # 790 Pretoria
1.7459557262705000E-15, # 372 Palma
1.6717209917006440E-15, # 107 Camilla
1.5850986571596890E-15, # 354 Eleonora
1.5465676956243250E-15, # 96 Aegle
1.5079333711965190E-15, # 386 Siegena
1.4975196825567010E-15, # 39 Laetitia
1.3886265898561990E-15, # 324 Bamberga
1.3315362554599750E-15, # 11 Parthenope
1.2792300000000000E-15, # 94 Aurora
1.2026244434834600E-15, # 22 Kalliope
1.1889849499520080E-15, # 120 Lachesis
1.1355858944839220E-15, # 185 Eunike
1.1006456795750680E-15, # 14 Irene
1.0975631032822510E-15, # 536 Merapi
1.0778410042407300E-15, # 9 Metis
1.0684599027046930E-15, # 165 Loreley
1.0356448401311940E-15, # 19 Fortuna
9.9366295459092490E-16, # 130 Elektra
9.6501295105487520E-16, # 128 Nemesis
9.4193644635953100E-16, # 187 Lamberta
9.3242237621988690E-16, # 41 Daphne
9.3159485940656200E-16, # 532 Herculina
9.2540854530185370E-16, # 85 Io
9.0708048441145050E-16, # 444 Gyptis
8.8950672849270450E-16, # 702 Alauda
8.5612605995538920E-16, # 137 Meliboea
8.4594307289596830E-16, # 471 Papagena
8.4256780185679340E-16, # 45 Eugenia
8.3124192126733720E-16, # 6 Hebe
8.3122000000000010E-16, # 76 Freia
8.2983887671636940E-16, # 154 Bertha
8.2144999999999990E-16, # 409 Aspasia
7.5582329262288920E-16, # 145 Adeona
7.2405478852581390E-16, # 308 Polyxo
7.0079069220413430E-16, # 121 Hermione
7.0078739271302920E-16, # 349 Dembowska
6.9951633549830870E-16, # 144 Vibilia
6.9079712474674250E-16, # 216 Kleopatra
6.3394427275876510E-16, # 59 Elpis
6.2808585539363830E-16, # 259 Aletheia
6.2674000000000000E-16, # 566 Stereoskopia
6.2157460662366880E-16, # 747 Winchester
5.8942565297069080E-16, # 8 Flora
5.6477307179764760E-16, # 93 Minerva
5.6241736501924590E-16, # 54 Alexandra
5.5764804768085390E-16, # 405 Thia
5.5435235615988890E-16, # 47 Aglaja
5.4846724291131710E-16, # 489 Comacina
5.3973999999999980E-16, # 106 Dione
5.3689097042583350E-16, # 344 Desiderata
5.3470999999999980E-16, # 168 Sibylla
5.2966693580739110E-16, # 238 Hypatia
5.2669011092543480E-16, # 268 Adorea
5.1446100208767350E-16, # 69 Hesperia
5.0204242212925650E-16, # 712 Boliviana
4.8352000000000000E-16, # 420 Bertholda
4.8122396678018730E-16, # 104 Klymene
4.7637535316542470E-16, # 690 Wratislavia
4.6542452473979750E-16, # 129 Antigone
4.6386880160793900E-16, # 618 Elfriede
4.5585992488069250E-16, # 375 Ursula
4.5559595663772370E-16, # 814 Tauris
4.5336468487968120E-16, # 150 Nuwa
4.5014136561179280E-16, # 196 Philomela
4.4713680178417890E-16, # 117 Lomia
4.2242882143774450E-16, # 139 Juewa
4.1832602031724900E-16, # 164 Eva
4.1456571725352780E-16, # 225 Henrietta
4.0681729177164310E-16, # 147 Protogeneia
4.0369435176860730E-16, # 92 Undina
3.9416000000000000E-16, # 146 Lucina
3.9124285470831010E-16, # 173 Ino
3.8880038985457880E-16, # 27 Euterpe
3.8674996577516790E-16, # 212 Medea
3.8350891708780020E-16, # 596 Scheila
3.7864750330894840E-16, # 895 Helio
3.7661406016708330E-16, # 141 Lumen
3.7487362845520320E-16, # 5 Astraea
3.7154667975341450E-16, # 105 Artemis
3.6992883127021260E-16, # 56 Melete
3.6806019206396510E-16, # 57 Mnemosyne
3.6782000000000000E-16, # 419 Aurelia
3.6611682430617230E-16, # 127 Johanna
3.6414973977832910E-16, # 490 Veritas
3.6186546970297340E-16, # 410 Chloris
3.5953000000000000E-16, # 654 Zelinda
3.5934934807617610E-16, # 381 Myrrha
3.5073744512956140E-16, # 74 Galatea
3.4333444128591690E-16, # 388 Charybdis
3.4310265912379690E-16, # 68 Leto
3.4176985017472390E-16, # 505 Cava
3.4106769797534870E-16, # 508 Princetonia
3.4020311574394290E-16, # 89 Julia
3.3677605743530070E-16, # 360 Carlova
3.3629754927610390E-16, # 909 Ulla
3.3620465420887160E-16, # 134 Sophrosyne
3.3407005297045190E-16, # 481 Emita
3.2728000000000000E-16, # 46 Hestia
3.2624933794345620E-16, # 334 Chicago
3.2284000000000000E-16, # 469 Argentina
3.1316732532224060E-16, # 140 Siwa
3.0567119446636530E-16, # 776 Berbericia
3.0465115564904100E-16, # 211 Isolda
3.0054836259460050E-16, # 241 Germania
2.9524140803084220E-16, # 40 Harmonia
2.9445412915212860E-16, # 34 Circe
2.9376520795314820E-16, # 514 Armida
2.9262727442945280E-16, # 28 Bellona
2.9175003506482500E-16, # 705 Erminia
2.9127000000000000E-16, # 328 Gudrun
2.8990971962115820E-16, # 175 Andromache
2.8848684039510110E-16, # 303 Josephina
2.8824066183805000E-16, # 159 Aemilia
2.8716404826701960E-16, # 1093 Freda
2.7688888401578460E-16, # 70 Panopaea
2.7653166434743810E-16, # 42 Isis
2.7485668803401500E-16, # 554 Peraga
2.7230587289059820E-16, # 194 Prokne
2.7166197083932590E-16, # 95 Arethusa
2.7097270724162130E-16, # 247 Eukrate
2.7071416736527810E-16, # 466 Tisiphone
2.6816131961351820E-16, # 356 Liguria
2.6336000000000000E-16, # 156 Xanthippe
2.5931073098872930E-16, # 209 Dido
2.5705491133531450E-16, # 51 Nemausa
2.5294428720409990E-16, # 21 Lutetia
2.5206144925338950E-16, # 772 Tanete
2.5119412656578330E-16, # 192 Nausikaa
2.4428317417320690E-16, # 98 Ianthe
2.4404616777010060E-16, # 91 Aegina
2.4193160564640660E-16, # 476 Hedwig
2.3729075379342010E-16, # 162 Laurentia
2.3522561732418410E-16, # 35 Leukothea
2.2979758127971450E-16, # 426 Hippo
2.2861593339578100E-16, # 455 Bruchsalia
2.2671033414599050E-16, # 804 Hispania
2.2481550488436830E-16, # 276 Adelheid
2.2211133615828450E-16, # 595 Polyxena
2.2099166077176560E-16, # 346 Hermentaria
2.1856205771130560E-16, # 37 Fides
2.1581556817636010E-16, # 788 Hohensteina
2.1523995570228910E-16, # 86 Semele
2.1501467433612250E-16, # 602 Marianna
2.0818012830905660E-16, # 762 Pulcova
2.0815063964697380E-16, # 17 Thetis
2.0528101779869580E-16, # 283 Emma
2.0089277366511320E-16, # 18 Melpomene
1.9735313911048410E-16, # 506 Marion
1.9725153985693760E-16, # 769 Tatjana
1.9715919666254550E-16, # 233 Asterope
1.9688550185992420E-16, # 250 Bettina
1.9576160906423330E-16, # 773 Irmintraud
1.9484986115571250E-16, # 545 Messalina
1.9317757851829200E-16, # 12 Victoria
1.9151562798850780E-16, # 488 Kreusa
1.8953317604197830E-16, # 23 Thalia
1.8915247466562090E-16, # 326 Tamara
1.8654000000000000E-16, # 275 Sapientia
1.8492472080300540E-16, # 203 Pompeja
1.8413301875237820E-16, # 266 Aline
1.8168406201733030E-16, # 221 Eos
1.8124493964458610E-16, # 521 Brixia
1.7989374811476500E-16, # 751 Faina
1.7619560261002570E-16, # 357 Ninina
1.7608987155151350E-16, # 230 Athamantis
1.7554616782526140E-16, # 200 Dynamene
1.7252555236382890E-16, # 191 Kolga
1.7050000000000000E-16, # 114 Kassandra
1.6995669526331390E-16, # 635 Vundtia
1.6970600184097090E-16, # 36 Atalante
1.6796766833548250E-16, # 181 Eucharis
1.6713475520093610E-16, # 176 Iduna
1.6643778725882150E-16, # 626 Notburga
1.6595259516345300E-16, # 148 Gallia
1.6373439522610840E-16, # 26 Proserpina
1.6333263911175180E-16, # 50 Virginia
1.6308230966598090E-16, # 784 Pickeringia
1.5923470884492270E-16, # 675 Ludmilla
1.5830198804699920E-16, # 786 Bredichina
1.5678504140474830E-16, # 407 Arachne
1.5586000000000000E-16, # 393 Lampetia
1.5456008538659910E-16, # 709 Fringilla
1.5352978885566200E-16, # 980 Anacostia
1.5109884909293550E-16, # 171 Ophelia
1.4938668030133560E-16, # 694 Ekard
1.4869856293449610E-16, # 416 Vaticana
1.4858970838252890E-16, # 491 Carina
1.4857405347656310E-16, # 780 Armenia
1.4820190164375290E-16, # 30 Urania
1.4756769258889110E-16, # 335 Roberta
1.4655734975631910E-16, # 412 Elisabetha
1.4521678954842190E-16, # 404 Arsinoe
1.4383748004467890E-16, # 674 Rachele
1.4268157848333010E-16, # 713 Luscinia
1.4244927463509560E-16, # 71 Niobe
1.4144228920898080E-16, # 377 Campania
1.4100078405763090E-16, # 350 Ornamenta
1.4076985722105040E-16, # 110 Lydia
1.3906039685198460E-16, # 223 Rosa
1.3720617596137790E-16, # 373 Melusina
1.3659771964686970E-16, # 100 Hekate
1.3558171306734880E-16, # 449 Hamburga
1.3232859647467680E-16, # 38 Leda
1.3145378148802100E-16, # 210 Isabella
1.2987929239170220E-16, # 498 Tokio
1.2812662256605980E-16, # 102 Miriam
1.2583765926450050E-16, # 1015 Christa
1.2573126556318860E-16, # 84 Klio
1.2351963628284910E-16, # 90 Antiope
1.2317521169884130E-16, # 345 Tercidina
1.2196487517410750E-16, # 791 Ani
1.2178097995831790E-16, # 358 Apollonia
1.2049549488578180E-16, # 663 Gerlinde
1.2048459545465020E-16, # 227 Philosophia
1.1995850162334400E-16, # 32 Pomona
1.1892814439317690E-16, # 740 Cantabia
1.1789669849357310E-16, # 206 Hersilia
1.1614439541131080E-16, # 80 Sappho
1.1584591337407960E-16, # 313 Chaldaea
1.1563954713645830E-16, # 366 Vincentina
1.1534210849315850E-16, # 445 Edna
1.1363293901133810E-16, # 236 Honoria
1.1260861147488830E-16, # 143 Adria
1.1187183314500760E-16, # 503 Evelyn
1.1177734423974850E-16, # 696 Leonora
1.1152801330348170E-16, # 1467 Mashona
1.1048622528732650E-16, # 240 Vanadis
1.1031851822968280E-16, # 385 Ilmatar
1.0968348906260350E-16, # 83 Beatrix
1.0890481919600570E-16, # 62 Erato
1.0886266508800200E-16, # 517 Edith
1.0826185861581930E-16, # 109 Felicitas
1.0436103612506560E-16, # 683 Lanzia
1.0354972501527590E-16, # 849 Ara
1.0314956358376310E-16, # 97 Klotho
1.0248487058744050E-16, # 160 Una
1.0223675545561340E-16, # 81 Terpsichore
1.0147296895716560E-16, # 201 Penelope
1.0044659839630940E-16, # 387 Aquitania
1.0011180168586460E-16, # 103 Hera
9.9581620100255190E-17, # 680 Genoveva
9.5491282158876470E-17, # 1107 Lictoria
9.5223891872172390E-17, # 213 Lilaea
9.5156050418462080E-17, # 135 Hertha
9.4537052258201020E-17, # 691 Lehigh
9.1998074776309120E-17, # 20 Massalia
9.1514999225958940E-17, # 205 Martha
9.0942944157907600E-17, # 667 Denise
8.8772708265623380E-17, # 124 Alkeste
8.6580000000000000E-17, # 163 Erigone
8.5978021881831310E-17, # 568 Cheruskia
8.5976202637173240E-17, # 735 Marghanna
8.4811739114660020E-17, # 58 Concordia
8.4019062534638870E-17, # 78 Diana
8.3548198507378080E-17, # 195 Eurykleia
8.3518243314079400E-17, # 79 Eurynome
8.2763495213043580E-17, # 322 Phaeo
8.0870792881032830E-17, # 362 Havnia
8.0635184290055780E-17, # 464 Megaira
7.9950510449165410E-17, # 72 Feronia
7.7695552096218000E-17, # 365 Corduba
7.6808847746999990E-17, # 454 Mathesis
7.6171409872842650E-17, # 1021 Flammario
7.5494816293144010E-17, # 49 Pales
7.3799853443752290E-17, # 424 Gratia
7.3526620138415910E-17, # 99 Dike
7.2753938533407130E-17, # 43 Ariadne
7.2538616434661140E-17, # 535 Montague
7.2398415223662110E-17, # 25 Phocaea
7.1382433975293480E-17, # 739 Mandeville
6.9841945302118660E-17, # 363 Padua
6.9600000000000000E-17, # 516 Amherstia
6.9118774732740100E-17, # 599 Luisa
6.8832710966054010E-17, # 304 Olga
6.6012607669307700E-17, # 82 Alkmene
6.5399257724402350E-17, # 465 Alekto
6.4677672132374610E-17, # 1036 Ganymed
6.4260399533186430E-17, # 604 Tekmessa
6.3407783495216970E-17, # 431 Nephele
6.2485687084637610E-17, # 1171 Rusthawelia
6.2392433107751650E-17, # 53 Kalypso
6.2251018387380330E-17, # 389 Industria
6.1530641540755240E-17, # 598 Octavia
6.1123219726683140E-17, # 569 Misa
5.8167921005732780E-17, # 760 Massinga
5.7960397015532350E-17, # 112 Iphigenia
5.7345815663368580E-17, # 369 Aeria
5.6430000000000000E-17, # 336 Lacadiera
5.6404017439762430E-17, # 63 Ausonia
5.5534656287469120E-17, # 442 Eichsfeldia
5.5258241903855270E-17, # 115 Thyra
5.4001214087434230E-17, # 415 Palatia
5.1407683428632200E-17, # 338 Budrosa
5.1143093942960570E-17, # 329 Svea
5.0911367830144640E-17, # 60 Echo
5.0188993452074920E-17, # 593 Titania
4.9982048125892480E-17, # 752 Sulamitis
4.9312955095007290E-17, # 77 Frigga
4.8033737856096100E-17, # 287 Nephthys
4.7743547738799830E-17, # 778 Theobalda
4.6886407201292200E-17, # 44 Nysa
4.5372886294250990E-17, # 177 Irma
4.3573746250771270E-17, # 75 Eurydike
4.3177157644920530E-17, # 172 Baucis
4.1172682683995420E-17, # 914 Palisana
3.9520590343452720E-17, # 224 Oceana
3.6737960799234180E-17, # 485 Genua
3.4515074511865960E-17, # 337 Devosa
3.3519192811280560E-17, # 111 Ate
2.8123577458657680E-17, # 547 Praxedis
2.7000759625981350E-17, # 118 Peitho
2.5580213924578190E-17, # 113 Amalthea
2.1949122177817370E-17, # 347 Pariana
2.1760595614159460E-17, # 584 Semiramis
2.0538576829083310E-17, # 198 Ampella
2.0038393509521630E-17, # 591 Irmgard
1.8780122958989280E-17, # 432 Pythia
1.4497976797693290E-17, # 623 Chimaera
1.3196141221701560E-17, # 132 Aethra
1.1716991540679400E-17, # 585 Bilkis
9.9000011897959020E-19  # 433 Eros
]

# Vector of JPL DE430 asteroids IDs
# See Table 12 in pages 53-59 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const de430_343ast_ids =
    [1, 4, 2, 10, 31, 704, 511, 15, 3, 16,
    65, 88, 48, 52, 451, 87, 7, 423, 29, 24,
    13, 790, 372, 107, 354, 96, 386, 39, 324, 11,
    94, 22, 120, 185, 14, 536, 9, 165, 19, 130,
    128, 187, 41, 532, 85, 444, 702, 137, 471, 45,
    6, 76, 154, 409, 145, 308, 121, 349, 144, 216,
    59, 259, 566, 747, 8, 93, 54, 405, 47, 489,
    106, 344, 168, 238, 268, 69, 712, 420, 104, 690,
    129, 618, 375, 814, 150, 196, 117, 139, 164, 225,
    147, 92, 146, 173, 27, 212, 596, 895, 141, 5,
    105, 56, 57, 419, 127, 490, 410, 654, 381, 74,
    388, 68, 505, 508, 89, 360, 909, 134, 481, 46,
    334, 469, 140, 776, 211, 241, 40, 34, 514, 28,
    705, 328, 175, 303, 159, 1093, 70, 42, 554, 194,
    95, 247, 466, 356, 156, 209, 51, 21, 772, 192,
    98, 91, 476, 162, 35, 426, 455, 804, 276, 595,
    346, 37, 788, 86, 602, 762, 17, 283, 18, 506,
    769, 233, 250, 773, 545, 12, 488, 23, 326, 275,
    203, 266, 221, 521, 751, 357, 230, 200, 191, 114,
    635, 36, 181, 176, 626, 148, 26, 50, 784, 675,
    786, 407, 393, 709, 980, 171, 694, 416, 491, 780,
    30, 335, 412, 404, 674, 713, 71, 377, 350, 110,
    223, 373, 100, 449, 38, 210, 498, 102, 1015, 84,
    90, 345, 791, 358, 663, 227, 32, 740, 206, 80,
    313, 366, 445, 236, 143, 503, 696, 1467, 240, 385,
    83, 62, 517, 109, 683, 849, 97, 160, 81, 201,
    387, 103, 680, 1107, 213, 135, 691, 20, 205, 667,
    124, 163, 568, 735, 58, 78, 195, 79, 322, 362,
    464, 72, 365, 454, 1021, 49, 424, 99, 43, 535,
    25, 739, 363, 516, 599, 304, 82, 465, 1036, 604,
    431, 1171, 53, 389, 598, 569, 760, 112, 369, 336,
    63, 442, 115, 415, 338, 329, 60, 593, 752, 77,
    287, 778, 44, 177, 75, 172, 914, 224, 485, 337,
    111, 547, 118, 113, 347, 584, 198, 591, 432, 623,
    132, 585, 433
]

# Matrix of second-degree zonal harmonic J_2 interactions included in DE430 ephemeris
# See Section III. Translational Equations of Motion in pages 11-15 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const UJ_interaction = fill(false, length(μ), length(μ))
# First paragraph of Section III. Translational Equations of Motion in page 11 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract:
# Sun's J_2 only interacts with the Moon and planets
UJ_interaction[2:11, su] .= true
# Earth's grav potential interacts with Sun, Mercury, Venus, Moon, Mars and Jupiter
UJ_interaction[union(1:ea-1,ea+1:7), ea] .= true 
# Moon's grav potential interacts with Sun, Mercury, Venus, Earth, Mars and Jupiter
UJ_interaction[union(1:mo-1,mo+1:7), mo] .= true 

# Physical constants in various units 
# See Table 4 in page 47 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const au = 1.495978707E8                       # Astronomical unit value in km
const yr = 365.25                              # Days in a Julian year
const daysec = 86_400                          # Number of seconds in a day

# Parameters related to speed of light, c
const clightkms = 2.99792458E5                 # Speed of light, km/sec
const c_au_per_day = daysec*(clightkms/au)     # Speed of light in au per day
const c_au_per_sec = clightkms/au              # Speed of light in au per sec
const c_cm_per_sec = 100_000*clightkms         # Speed of light in cm per sec
const c_p2 = 29979.063823897606                # Speed of light^2 in au^2/day^2
const c_m2 = 3.3356611996764786e-5             # Speed of light^-2 in day^2/au^2

const sundofs = nbodyind(length(μ), su)        # Sun's position and velocity indices
const earthdofs = nbodyind(length(μ), ea)      # Earth's position and velocity indices

const J2000 = 2451545.0                        # 2000 epoch julian date   

# Extended body parameters for the Sun.
# See Table 9 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const R_sun = 696000.0/au                      # Solar radius in au

const α_p_sun = 286.13                         # Sun's rotation pole right ascension (degrees)
const δ_p_sun = 63.87                          # Sun's rotation pole declination (degrees)

const RSUN  = 6.9600000000000000E+05           # Solar radius in km  
const J2SUN = 2.1106088532726840E-07           # Dynamical form–factor of the Sun
# Second zonal harmonic coefficient J_2 * Suns's radius in au^2
const JS = [J2SUN*(RSUN/au)^2]                 

# Extended body parameters for the Earth
# See Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const RE     =  6.378136300E+03                # Earth's radius in km 
# TODO: solve differences between parsed and non-parsed
const RE_au = (RE/au)                          # Earth's radius in au
const J2E    =  1.082625450E-03                # Second zonal harmonic of the Earth 
const J3E    = -2.532410000E-06                # Third zonal harmonic of the Earth 
const J4E    = -1.619898000E-06                # Fourth zonal harmonic of the Earth 
const J5E    = -2.277345000E-07                # Fifth zonal harmonic of the Earth 
const J2EDOT = -2.600000000E-11                # Rate of change of J2E in 1/yr
# Vector of zonal harmonic coefficients J_n * Earth's radius in au ^n 
const JE = [J2E*(RE/au)^2, J3E*(RE/au)^3, J4E*(RE/au)^4, J5E*(RE/au)^5]

# Extended body parameters for the Moon

# Extended body parameters for the Moon
# See Table 11 in page 51 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const RM = 1.7380000000000000E+03            # Lunar radius in km 
const R_moon = RM/au                         # Lunar radius in au
const β_L = 6.3102131934887270E-04           # Lunar moment parameter, β_L = (C_T-A_T)/B_T
const γ_L = 2.2773171480091860E-04           # Lunar moment parameter, γ_L = (B_T-A_T)/C_T
const k_2M = 0.024059                        # Potential Love number
const τ_M = 9.5830547273306690E-02           # Time-lag for the lunar solid-body tide (days)
const α_c = 0.0007                           # Ratio of polar moment of inertia of core to mean total polar moment of inertia
const f_c = 2.4623904789198150E-04           # Oblateness of core
const k_ν_div_C_T = 1.6365616533709530E-08   # Friction coefficient between core and mantle, radian/day

const J2_M_und = 2.0321568464952570E-04      # Undistorted lunar 2nd zonal harmonic coefficient

const J2M  =  2.0321568464952570E-04           # Undistorted 2nd zonal harmonic coefficient
const J3M  =  8.4597026974594570E-06           # Third zonal harmonic coefficient
const J4M  = -9.7044138365700000E-06           # Fourth zonal harmonic coefficient
const J5M  =  7.4221608384052890E-07           # Fifth zonal harmonic coefficient
const J6M  = -1.3767531350969900E-05           # Sixth zonal harmonic coefficient
# Vector of zonal harmonic coefficients J_n * Moons's radius in au ^n
const JM = [J2M*R_moon^2, J3M*R_moon^3, J4M*R_moon^4, J5M*R_moon^5, J6M*R_moon^6]

# Matrix of zonal harmonic coefficients J_n * radius of the corresponding body ^n
const JSEM = zeros(5, 6)
JSEM[su,2:2] = JS     # Sun 
JSEM[ea,2:5] = JE     # Earth 
JSEM[mo,2:6] = JM     # Moon

# Degree of zonal harmonics expansion for each body with extended body accelerations
const n1SEM = [2, 0, 0, 5, 6]   # Sun, Mercury, Venus, Earth, Moon
# Degree of zonal harmonics expansion for the Moon
const n2M = 6

# Numerical factors for recursion relations of Legendre polynomials
# See equations (175)-(178) in page 33 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
const max_n1SEM = maximum(n1SEM)
const fact1_jsem = [(2n-1)/n for n in 1:max_n1SEM]  # (2n - 1) / n
const fact2_jsem = [(n-1)/n for n in 1:max_n1SEM]   # (n - 1) / n
const fact3_jsem = [n for n in 1:max_n1SEM]         # n
const fact4_jsem = fact3_jsem .+ 1                  # n + 1
const fact5_jsem = fact3_jsem .+ 2                  # n + 2    

# Lunar tesseral harmonics coefficients (C_{nm}, S_{nm}) with n = 2,...,6 and m = 1,...,n
# See Table 11 in page 51 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract

# n = 2, m = 2
const C22M =  2.2382740590560020E-05          

# n = 3, m = 1, 2, 3
const C31M =  2.8480741195592860E-05
const S31M =  5.8915551555318640E-06
const C32M =  4.8449420619770600E-06
const S32M =  1.6844743962783900E-06
const C33M =  1.6756178134114570E-06
const S33M = -2.4742714379805760E-07
# Vector of tesseral harmonic coefficients (C_{3m}, S_{3m}) * Moons's radius in au ^3
const C3M = [C31M, C32M, C33M]*R_moon^3
const S3M = [S31M, S32M, S33M]*R_moon^3

# n = 4, m = 1, 2, 3, 4
const C41M = -5.7048697319733210E-06
const S41M =  1.5789202789245720E-06
const C42M = -1.5912271792977430E-06
const S42M = -1.5153915796731720E-06
const C43M = -8.0678881596778210E-08
const S43M = -8.0349266627431070E-07
const C44M = -1.2692158612216040E-07
const S44M =  8.2964257754075220E-08
# Vector of tesseral harmonic coefficients (C_{4m}, S_{4m}) * Moons's radius in au ^4
const C4M = [C41M, C42M, C43M, C44M]*R_moon^4
const S4M = [S41M, S42M, S43M, S44M]*R_moon^4

# n = 5, m = 1, 2, 3, 4, 5
const C51M = -8.6629769308983560E-07
const S51M = -3.5272289393243820E-06
const C52M =  7.1199537967353330E-07
const S52M =  1.7107886673430380E-07
const C53M =  1.5399750424904520E-08
const S53M =  2.8736257616334340E-07
const C54M =  2.1444704319218450E-08
const S54M =  5.2652110720146800E-10
const C55M =  7.6596153884006140E-09
const S55M = -6.7824035473995330E-09
# Vector of tesseral harmonic coefficients (C_{5m}, S_{5m}) * Moons's radius in au ^5
const C5M = [C51M, C52M, C53M, C54M, C55M]*R_moon^5
const S5M = [S51M, S52M, S53M, S54M, S55M]*R_moon^5

# n = 6, m = 1, 2, 3, 4, 5, 6
const C61M =  1.2024363601545920E-06
const S61M = -2.0453507141252220E-06
const C62M = -5.4703897324156850E-07
const S62M = -2.6966834353574270E-07
const C63M = -6.8785612757292010E-08
const S63M = -7.1063745295915780E-08
const C64M =  1.2915580402925160E-09
const S64M = -1.5361616966632300E-08
const C65M =  1.1737698784460500E-09
const S65M = -8.3465073195142520E-09
const C66M = -1.0913395178881540E-09
const S66M =  1.6844213702632920E-09
# Vector of tesseral harmonic coefficients (C_{6m}, S_{6m}) * Moons's radius in au ^6
const C6M = [C61M, C62M, C63M, C64M, C65M, C66M]*R_moon^6
const S6M = [S61M, S62M, S63M, S64M, S65M, S66M]*R_moon^6

# Matrix of tesseral harmonic coefficients C_{nm} * radius of the moon ^n
const CM = zeros(6, 6)
CM[3,1:3] .= C3M    # n = 3
CM[4,1:4] .= C4M    # n = 4
CM[5,1:5] .= C5M    # n = 5
CM[6,1:6] .= C6M    # n = 6
# Matrix of tesseral harmonic coefficients S_{nm} * radius of the moon ^n
const SM = zeros(6, 6)
SM[3,1:3] .= S3M    # n = 3
SM[4,1:4] .= S4M    # n = 4
SM[5,1:5] .= S5M    # n = 5
SM[6,1:6] .= S6M    # n = 6

# Numerical factors for recursion relations of Associated Legendre polynomials
# See equations (180)-(183) in pages 33-34 of https://ui.adsabs.harvard.edu/abs/1971mfdo.book.....M/abstract
const lnm1 = [(2n-1)/(n-m) for n in 1:6, m in 1:6]       # (2n - 1) / (n - m)
const lnm2 = [-(n+m-1)/(n-m) for n in 1:6, m in 1:6]     # -(n + m - 1) / (n - m)
const lnm3 = [-n for n in 1:6]                           # -n 
const lnm4 = [n+m for n in 1:6, m in 1:6]                # (n + m)
const lnm5 = [2n-1 for n in 1:6]                         # (2n - 1)
const lnm6 = [-(n+1) for n in 1:6]                       # -(n + 1)
const lnm7 = [m for m in 1:6]                            # m  

# Diagonal elements of undistorted lunar mantle moment of inertia
# See equation (37) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const A_T = (2(1-β_L*γ_L)/(2β_L-γ_L+β_L*γ_L))*J2_M_und
const B_T = (2(1+γ_L)/(2β_L-γ_L+β_L*γ_L))*J2_M_und
const C_T = (2(1+β_L)/(2β_L-γ_L+β_L*γ_L))*J2_M_und

const k_ν = k_ν_div_C_T*(C_T*μ[mo]*(R_moon^2))

# Diagonal elements of lunar core moment of inertia
# See equation (39) in page 17 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const A_c = α_c*C_T*(1-f_c)
const B_c = α_c*C_T*(1-f_c)
const C_c = α_c*C_T
# Numerical factor for the torque on the mantle due to the interaction between the core and mantle
# See equation (45) in page 18 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const C_c_m_A_c = (C_c - A_c)*(μ[mo]*R_moon^2)

# Lunar undistorted total moment of inertia
# See equation (36) in page 16 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const ITM_und = diagm([A_T, B_T, C_T]*μ[mo]*R_moon^2)

# Lunar core moment of inertia
# See equation (39) in page 17 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const I_c = diagm([A_c, B_c, C_c]*μ[mo]*R_moon^2)

# Lunar distance (km)
const ld = 384402.0 
# 384399.014 km  mean semi-major axis value was retrieved from:
# Table 3, page 5, J. G. Williams, D. H. Boggs, and W. M. Folkner (2013). 
# “DE430 Lunar Orbit, Physical Librations and Surface Coordinates,” JPL IOM 335-JW,
# DB,WF-20130722-016 (internal document)
# https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430_moon_coord.pdf
const n_moon = sqrt((μ[ea]+μ[mo])/((384399.014/au)^3)) # Lunar mean motion (rad/day)

# Love numbers for Earth
# See Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const k_20E = 0.335           # Potential Love number for long-period deformation
const k_21E = 0.32            # Potential Love number for diurnal deformation
const k_22E = 0.32            # Potential Love number for semi-diurnal deformation

# Orbital time-lags  for Earth
# See Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const τ_0p = 0.0640           # Orbital time-lag for long-period deformation, days
const τ_1p = -0.044           # Orbital time-lag for diurnal deformation, days --> DE430 report states -0.0044, but tech comments and Williams et al. (2006) state -0.044
const τ_2p = -0.1000          # Orbital time-lag for semi-diurnal deformation, days

# Rotational time-lags for Earth
# See Table 10 in page 50 of https://ui.adsabs.harvard.edu/abs/2014IPNPR.196C...1F%2F/abstract
const τ_0 = 0.0                    # Rotational time-lag for long-period deformation, days
const τ_1 = 7.3632190228041890E-03 # Rotational time-lag for diurnal deformation, days
const τ_2 = 2.5352978633388720E-03 # Rotational time-lag for semi-diurnal deformation, days

# Standard value of nominal mean angular velocity of Earth (rad/day), 
# See Explanatory Supplement to the Astronomical Almanac 2014, section 7.4.3.3,
# page 296: 7.2921151467e-5 rad/second
const ω_E = daysec*7.2921151467e-5 # (7.2921151467e-5 rad/sec)*daysec -> rad/day

# Earth/Moon mass ratio
const EMRAT = 8.1300569074190620E+01
