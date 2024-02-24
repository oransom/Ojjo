
shittopo=readtable("D:\ORCaS\Bellefield2\topo\Ojjo_Bellefield_2_20240108_231215-xml_CA-5_5-usft.csv");
shitflood=readtable("D:\ORCaS\Bellefield2\2023-12-21_Bellefield 2_CombinedWSETrimmed_XYZ.csv");

NWX=6532339.8488;
NWY=2190724.1751;
NEX=6539216.5608;
NEY=2190724.1751;
SWX=6532339.8488;
SWY=2186639.6832;
SEX=6539216.5608;	
SEY=2186639.6832;

shittopo.Var1(shittopo.Var1<SWX | shittopo.Var1>NEX)=NaN;
shittopo.Var2(shittopo.Var2<SEY | shittopo.Var2>NEY)=NaN;
shitflood.Var1(shitflood.Var1<SWX | shitflood.Var1>NEX)=NaN;
shitflood.Var2(shitflood.Var2<SEY | shitflood.Var2>NEY)=NaN;

shittopo=rmmissing(shittopo);
shitflood=rmmissing(shitflood);

writetable(shittopo,'topob_1-14.csv');
writetable(shitflood,'floodb_1-14.csv');