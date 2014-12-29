/*
 * apeTreeBuilder.cpp
 * Makes a tree object according to the commonly accepted ape phylogeny.
 * The data was taken from Perelman, Polina, et al. "A molecular phylogeny of living primates." PLoS genetics 7.3 (2011): e1001342.
 * Author: Evan Albright, Jack Hessel, Nao Hiranuma, and Cody Wang
 * Last modified: 2/16/2014
 */

#include "apeTreeBuilder.h"

Tree ApeTreeBuilder::makeTree(){
  Tree tree;
  edge_descriptor e;

  Species* lophocebus = new Species("LophocebusAterrimus");
  vertex_descriptor vlophocebus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vlophocebus, lophocebus);

  Species* dmy = new Species("dmy17");
  vertex_descriptor n17 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n17, dmy);
  
  e = add_edge(vlophocebus, n17, tree.g).first;
  put(edge_weight_t(), tree.g, e, 3.24);

  Species* papio = new Species("PapioPapio");
  vertex_descriptor vpapio = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vpapio, papio);
  
  e = add_edge(vpapio, n17, tree.g).first;
  put(edge_weight_t(), tree.g, e, 3.24);

  dmy = new Species("dmy19");
  vertex_descriptor n19 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n19, dmy);
  
  e = add_edge(n19, n17, tree.g).first;
  put(edge_weight_t(), tree.g, e, 0.82);

  Species* theropithecus = new Species("TheropithecusGelada");
  vertex_descriptor vtheropithecus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vtheropithecus, theropithecus);
  
  e = add_edge(n19, vtheropithecus, tree.g).first;
  put(edge_weight_t(), tree.g, e, 4.06);

  ///

  dmy = new Species("dmy14");
  vertex_descriptor n14 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n14, dmy);

  e = add_edge(n14, n19, tree.g).first;
  put(edge_weight_t(), tree.g, e, 2.61);

  Species* mandrillus = new Species("MandrillusSphinx");
  vertex_descriptor vmandrillus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vmandrillus, mandrillus);
  
  dmy = new Species("dmy23");
  vertex_descriptor n23 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n23, dmy);

  e = add_edge(vmandrillus, n23, tree.g).first;
  put(edge_weight_t(), tree.g, e, 4.85);

  Species* cercocebus = new Species("CercocebusChrysogaster");
  vertex_descriptor vcercocebus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vcercocebus, cercocebus);

  e = add_edge(vcercocebus, n23, tree.g).first;
  put(edge_weight_t(), tree.g, e, 4.85);

  e = add_edge(n14, n23, tree.g).first;
  put(edge_weight_t(), tree.g, e, 1.82);

  dmy = new Species("dmy24");
  vertex_descriptor n24 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n24, dmy);

  e = add_edge(n14, n24, tree.g).first;
  put(edge_weight_t(), tree.g, e, 1.46);

  Species* macaca = new Species("MacacaMulatta");
  vertex_descriptor vmacaca = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vmacaca, macaca);

  e = add_edge(vmacaca, n24, tree.g).first;
  put(edge_weight_t(), tree.g, e, 8.13);

  dmy = new Species("dmy41");
  vertex_descriptor n41 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n41, dmy);

  e = add_edge(n41, n24, tree.g).first;
  put(edge_weight_t(), tree.g, e, 3.37);

  /////cercopithecini
  Species* erythrocebus = new Species("ErythrocebusPatas");
  vertex_descriptor verythrocebus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, verythrocebus, erythrocebus);

  dmy = new Species("dmy39");
  vertex_descriptor n39 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n39, dmy);

  e = add_edge(n39, verythrocebus, tree.g).first;
  put(edge_weight_t(), tree.g, e, 4.95);

  Species* chlorocebus = new Species("ChlorocebusAethiops");
  vertex_descriptor vchlorocebus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vchlorocebus, chlorocebus);

  e = add_edge(n39, vchlorocebus, tree.g).first;
  put(edge_weight_t(), tree.g, e, 4.95);

  dmy = new Species("dmy25");
  vertex_descriptor n25 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n25, dmy);

  e = add_edge(n39, n25, tree.g).first;
  put(edge_weight_t(), tree.g, e, 3.27);
  
  Species* cercopithecus = new Species("CercopithecusAlbogularis");
  vertex_descriptor vcercopithecus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vcercopithecus, cercopithecus);

  e = add_edge(vcercopithecus, n25, tree.g).first;
  put(edge_weight_t(), tree.g, e, 8.22);

  e = add_edge(n41, n25, tree.g).first;
  put(edge_weight_t(), tree.g, e, 3.28);
 
  dmy = new Species("dmy62");
  vertex_descriptor n62 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n62, dmy);

  e = add_edge(n41, n62, tree.g).first;
  put(edge_weight_t(), tree.g, e, 17.57-11.50);

  ////colobinae
  Species* pygathrix = new Species("PygathrixCinerea");
  vertex_descriptor vpygathrix = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vpygathrix, pygathrix);

  dmy = new Species("dmy54");
  vertex_descriptor n54 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n54, dmy);

  e = add_edge(vpygathrix, n54, tree.g).first;
  put(edge_weight_t(), tree.g, e, 6.21);
 
  Species* nasalis = new Species("NasalisIarvatus");
  vertex_descriptor vnasalis = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vnasalis, nasalis);

  e = add_edge(vnasalis, n54, tree.g).first;
  put(edge_weight_t(), tree.g, e, 6.21);

  dmy = new Species("dmy55");
  vertex_descriptor n55 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n55, dmy);

  e = add_edge(n55, n54, tree.g).first;
  put(edge_weight_t(), tree.g, e, 6.69-6.21);
  
  Species* rhinopithecus = new Species("RhinopithecusRoxellana");
  vertex_descriptor vrhinopithecus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vrhinopithecus, rhinopithecus);

  e = add_edge(n55, vrhinopithecus, tree.g).first;
  put(edge_weight_t(), tree.g, e, 6.69);

  dmy = new Species("dmy57");
  vertex_descriptor n57 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n57, dmy);

  e = add_edge(n55, n57, tree.g).first;
  put(edge_weight_t(), tree.g, e, 8.30-6.69);
  
  Species* semnopithecus = new Species("SemnopithecusEntellus");
  vertex_descriptor vsemnopithecus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vsemnopithecus, semnopithecus);

  dmy = new Species("dmy51");
  vertex_descriptor n51 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n51, dmy);

  e = add_edge(n51, vsemnopithecus, tree.g).first;
  put(edge_weight_t(), tree.g, e, 4.05);
  
  Species* trachypithecus = new Species("TrachypithecusJohnii");
  vertex_descriptor vtrachypithecus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vtrachypithecus, trachypithecus);

  e = add_edge(n51, vtrachypithecus, tree.g).first;
  put(edge_weight_t(), tree.g, e, 4.05);

  e = add_edge(n51, n57, tree.g).first;
  put(edge_weight_t(), tree.g, e, 8.30-4.05);

  dmy = new Species("dmy58");
  vertex_descriptor n58 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n58, dmy);

  e = add_edge(n58, n57, tree.g).first;
  put(edge_weight_t(), tree.g, e, 8.81-8.30);
  
  Species* presbytis = new Species("PresbytisMelalophos");
  vertex_descriptor vpresbytis = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vpresbytis, presbytis);

  e = add_edge(n58, vpresbytis, tree.g).first;
  put(edge_weight_t(), tree.g, e, 8.81);

  dmy = new Species("dmy42");
  vertex_descriptor n42 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n42, dmy);

  e = add_edge(n58, n42, tree.g).first;
  put(edge_weight_t(), tree.g, e, 12.28-8.81);
  
  Species* piliocolobus = new Species("PiliocolobusBadius");
  vertex_descriptor vpiliocolobus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vpiliocolobus, piliocolobus);

  dmy = new Species("dmy61");
  vertex_descriptor n61 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n61, dmy);

  e = add_edge(n61, vpiliocolobus, tree.g).first;
  put(edge_weight_t(), tree.g, e, 9.12);
  
  Species* colobus = new Species("ColobusGuereza");
  vertex_descriptor vcolobus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vcolobus, colobus);

  e = add_edge(n61, vcolobus, tree.g).first;
  put(edge_weight_t(), tree.g, e, 9.12);

  dmy = new Species("dmy77");
  vertex_descriptor n77 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n77, dmy);

  e = add_edge(n61, n42, tree.g).first;
  put(edge_weight_t(), tree.g, e, 12.28-9.12);

  e = add_edge(n62, n42, tree.g).first;
  put(edge_weight_t(), tree.g, e, 17.57-12.28);

  e = add_edge(n62, n77, tree.g).first;
  put(edge_weight_t(), tree.g, e, 31.56-17.57);
  
  ////HOMONIDAE
  Species* pan = new Species("PanTroglodytes");
  vertex_descriptor vpan = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vpan, pan);

  dmy = new Species("dmy73");
  vertex_descriptor n73 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n73, dmy);

  e = add_edge(n73, vpan, tree.g).first;
  put(edge_weight_t(), tree.g, e, 6.60);
  
  Species* homo = new Species("HomoSapiens");
  vertex_descriptor vhomo = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vhomo, homo);

  e = add_edge(n73, vhomo, tree.g).first;
  put(edge_weight_t(), tree.g, e, 6.60);

  dmy = new Species("dmy74");
  vertex_descriptor n74 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n74, dmy);

  e = add_edge(n73, n74, tree.g).first;
  put(edge_weight_t(), tree.g, e, 8.30-6.60);
  
  Species* gorilla = new Species("GorillaGorilla");
  vertex_descriptor vgorilla = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vgorilla, gorilla);

  e = add_edge(vgorilla, n74, tree.g).first;
  put(edge_weight_t(), tree.g, e, 8.30);

  dmy = new Species("dmy63");
  vertex_descriptor n63 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n63, dmy);

  e = add_edge(n63, n74, tree.g).first;
  put(edge_weight_t(), tree.g, e, 20.32-8.30);
  
  Species* hylobates = new Species("HylobatesIar");
  vertex_descriptor vhylobates = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vhylobates, hylobates);

  dmy = new Species("dmy66");
  vertex_descriptor n66 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n66, dmy);

  e = add_edge(n66, vhylobates, tree.g).first;
  put(edge_weight_t(), tree.g, e, 8.00);
  
  Species* symphalangus = new Species("SymphalangusSyndactylus");
  vertex_descriptor vsymphalangus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vsymphalangus, symphalangus);

  e = add_edge(n66, vsymphalangus, tree.g).first;
  put(edge_weight_t(), tree.g, e, 8.00);

  dmy = new Species("dmy70");
  vertex_descriptor n70 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n70, dmy);

  e = add_edge(n66, n70, tree.g).first;
  put(edge_weight_t(), tree.g, e, 8.93-8.00);
  
  Species* nomascus = new Species("NomascusLeucogenys");
  vertex_descriptor vnomascus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vnomascus, nomascus);

  e = add_edge(vnomascus, n70, tree.g).first;
  put(edge_weight_t(), tree.g, e, 8.93);

  e = add_edge(n63, n70, tree.g).first;
  put(edge_weight_t(), tree.g, e, 20.32-8.93);

  e = add_edge(n63, n77, tree.g).first;
  put(edge_weight_t(), tree.g, e, 31.56-20.32);

  dmy = new Species("dmy141");
  vertex_descriptor n141 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n141, dmy);

  e = add_edge(n141, n77, tree.g).first;
  put(edge_weight_t(), tree.g, e, 42.47-31.56);

  //New world monkeys
  //Cebidae
  Species* callithrix = new Species("CallithrixPygmaea");
  vertex_descriptor vcallithrix = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vcallithrix, callithrix);

  dmy = new Species("dmy88");
  vertex_descriptor n88 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n88, dmy);

  e = add_edge(vcallithrix, n88, tree.g).first;
  put(edge_weight_t(), tree.g, e, 13.55);
  
  Species* leontopithecus = new Species("LeontopithecusRosalia");
  vertex_descriptor vleontopithecus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vleontopithecus, leontopithecus);

  e = add_edge(vleontopithecus, n88, tree.g).first;
  put(edge_weight_t(), tree.g, e, 13.55);

  dmy = new Species("dmy97");
  vertex_descriptor n97 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n97, dmy);

  e = add_edge(n97, n88, tree.g).first;
  put(edge_weight_t(), tree.g, e, 14.89-13.55);
  
  Species* saguinus = new Species("SaguinusOedipus");
  vertex_descriptor vsaguinus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vsaguinus, saguinus);

  e = add_edge(n97, vsaguinus, tree.g).first;
  put(edge_weight_t(), tree.g, e, 14.89);

  dmy = new Species("dmy112");
  vertex_descriptor n112 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n112, dmy);

  e = add_edge(n97, n112, tree.g).first;
  put(edge_weight_t(), tree.g, e, 19.25-14.89);
  
  Species* aotus = new Species("AotusNancymaae");
  vertex_descriptor vaotus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vaotus, aotus);

  e = add_edge(vaotus, n112, tree.g).first;
  put(edge_weight_t(), tree.g, e, 19.25);

  dmy = new Species("dmy113");
  vertex_descriptor n113 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n113, dmy);

  e = add_edge(n113, n112, tree.g).first;
  put(edge_weight_t(), tree.g, e, 19.95-19.25);
  
  Species* cebus = new Species("CebusXanthosternos");
  vertex_descriptor vcebus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vcebus, cebus);

  dmy = new Species("dmy111");
  vertex_descriptor n111 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n111, dmy);

  e = add_edge(n111, vcebus, tree.g).first;
  put(edge_weight_t(), tree.g, e, 15.40);
  
  Species* saimiri = new Species("SaimiriOerstedii");
  vertex_descriptor vsaimiri = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vsaimiri, saimiri);

  e = add_edge(n111, vsaimiri, tree.g).first;
  put(edge_weight_t(), tree.g, e, 15.40);

  e = add_edge(n113, n111, tree.g).first;
  put(edge_weight_t(), tree.g, e, 19.95-15.40);

  dmy = new Species("dmy127");
  vertex_descriptor n127 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n127, dmy);

  e = add_edge(n113, n127, tree.g).first;
  put(edge_weight_t(), tree.g, e, 22.76-19.95);
  
  Species* lagothrix = new Species("LagothrixLagotricha");
  vertex_descriptor vlagothrix = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vlagothrix, lagothrix);

  dmy = new Species("dmy122");
  vertex_descriptor n122 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n122, dmy);

  e = add_edge(n122, vlagothrix, tree.g).first;
  put(edge_weight_t(), tree.g, e, 11.25);
  
  Species* ateles = new Species("AtelesBelzebuth");
  vertex_descriptor vateles = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vateles, ateles);

  e = add_edge(n122, vateles, tree.g).first;
  put(edge_weight_t(), tree.g, e, 11.25);

  dmy = new Species("dmy126");
  vertex_descriptor n126 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n126, dmy);

  e = add_edge(n122, n126, tree.g).first;
  put(edge_weight_t(), tree.g, e, 16.13-11.25);
  
  Species* alouatta = new Species("AlouattaCaraya");
  vertex_descriptor valouatta = add_vertex(tree.g);
  put(Vertex_t(), tree.g, valouatta, alouatta);

  e = add_edge(valouatta, n126, tree.g).first;
  put(edge_weight_t(), tree.g, e, 16.13);

  dmy = new Species("dmy78");
  vertex_descriptor n78 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n78, dmy);

  e = add_edge(n127, n126, tree.g).first;
  put(edge_weight_t(), tree.g, e, 22.76-16.13);

  e = add_edge(n78, n127, tree.g).first;
  put(edge_weight_t(), tree.g, e, 24.82-22.76);

  e = add_edge(n141, n78, tree.g).first;
  put(edge_weight_t(), tree.g, e, 43.47-24.82);

  Species* cacajao = new Species("CacajaoCalvus");
  vertex_descriptor vcacajao = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vcacajao, cacajao);

  dmy = new Species("dmy137");
  vertex_descriptor n137 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n137, dmy);

  e = add_edge(n137, vcacajao, tree.g).first;
  put(edge_weight_t(), tree.g, e, 7.51);
  
  Species* chiropotes = new Species("ChiropotesAlbinasus");
  vertex_descriptor vchiropotes = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vchiropotes, chiropotes);

  e = add_edge(n137, vchiropotes, tree.g).first;
  put(edge_weight_t(), tree.g, e, 7.51);

  dmy = new Species("dmy140");
  vertex_descriptor n140 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n140, dmy);

  e = add_edge(n137, n140, tree.g).first;
  put(edge_weight_t(), tree.g, e, 20.24-7.51);
  
  Species* callicebus = new Species("CallicebusCupreus");
  vertex_descriptor vcallicebus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vcallicebus, callicebus);

  e = add_edge(vcallicebus, n140, tree.g).first;
  put(edge_weight_t(), tree.g, e, 20.24);

  e = add_edge(n78, n140, tree.g).first;
  put(edge_weight_t(), tree.g, e, 24.82-20.24);

  dmy = new Species("dmy143");
  vertex_descriptor n143 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n143, dmy);

  e = add_edge(n143, n141, tree.g).first;
  put(edge_weight_t(), tree.g, e, 81.00-43.47);
  
  Species* tarsius = new Species("TarsiusBancanus");
  vertex_descriptor vtarsius = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vtarsius, tarsius);

  e = add_edge(n143, vtarsius, tree.g).first;
  put(edge_weight_t(), tree.g, e, 81.00);

  dmy = new Species("dmy185");
  vertex_descriptor n185 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n185, dmy);

  e = add_edge(n143, n185, tree.g).first;
  put(edge_weight_t(), tree.g, e, 87.00-81.00);
  
  Species* lemur = new Species("LemurCatta");
  vertex_descriptor vlemur = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vlemur, lemur);

  dmy = new Species("dmy169");
  vertex_descriptor n169 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n169, dmy);

  e = add_edge(n169, vlemur, tree.g).first;
  put(edge_weight_t(), tree.g, e, 9.66);
  
  Species* hapalemur = new Species("HapalemurGriseus");
  vertex_descriptor vhapalemur = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vhapalemur, hapalemur);

  e = add_edge(n169, vhapalemur, tree.g).first;
  put(edge_weight_t(), tree.g, e, 9.66);

  dmy = new Species("dmy170");
  vertex_descriptor n170 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n170, dmy);

  e = add_edge(n169, n170, tree.g).first;
  put(edge_weight_t(), tree.g, e, 21.28-9.66);
  
  Species* eulemur = new Species("EulemurMacacoMacaco");
  vertex_descriptor veulemur = add_vertex(tree.g);
  put(Vertex_t(), tree.g, veulemur, eulemur);

  e = add_edge(veulemur, n170, tree.g).first;
  put(edge_weight_t(), tree.g, e, 21.28);

  dmy = new Species("dmy172");
  vertex_descriptor n172 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n172, dmy);

  e = add_edge(n172, n170, tree.g).first;
  put(edge_weight_t(), tree.g, e, 26.19-21.28);
  
  Species* varecia = new Species("VareciaVariegataVariegata");
  vertex_descriptor vvarecia = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vvarecia, varecia);

  e = add_edge(n172, vvarecia, tree.g).first;
  put(edge_weight_t(), tree.g, e, 26.19);

  dmy = new Species("dmy173");
  vertex_descriptor n173 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n173, dmy);

  e = add_edge(n172, n173, tree.g).first;
  put(edge_weight_t(), tree.g, e, 38.64-26.19);
  
  Species* cheirogaleus = new Species("CheirogaleusMedius");
  vertex_descriptor vcheirogaleus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vcheirogaleus, cheirogaleus);

  dmy = new Species("dmy152");
  vertex_descriptor n152 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n152, dmy);

  e = add_edge(vcheirogaleus, n152, tree.g).first;
  put(edge_weight_t(), tree.g, e, 32.92);
  
  Species* lepilemur = new Species("LepilemurHubbardorum");
  vertex_descriptor vlepilemur = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vlepilemur, lepilemur);

  e = add_edge(vlepilemur, n152, tree.g).first;
  put(edge_weight_t(), tree.g, e, 32.92);

  dmy = new Species("dmy158");
  vertex_descriptor n158 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n158, dmy);

  e = add_edge(n158, n152, tree.g).first;
  put(edge_weight_t(), tree.g, e, 37.44-32.92);
  
  Species* propithecus = new Species("PropithecusCoquereli");
  vertex_descriptor vpropithecus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vpropithecus, propithecus);

  dmy = new Species("dmy157");
  vertex_descriptor n157 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n157, dmy);

  e = add_edge(vpropithecus, n157, tree.g).first;
  put(edge_weight_t(), tree.g, e, 17.43);
  
  Species* avahi = new Species("AvahiLaniger");
  vertex_descriptor vavahi = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vavahi, avahi);

  e = add_edge(vavahi, n157, tree.g).first;
  put(edge_weight_t(), tree.g, e, 17.43);

  e = add_edge(n158, n157, tree.g).first;
  put(edge_weight_t(), tree.g, e, 37.44-17.43);

  e = add_edge(n158, n173, tree.g).first;
  put(edge_weight_t(), tree.g, e, 38.64-37.44);

  dmy = new Species("dmy174");
  vertex_descriptor n174 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n174, dmy);

  e = add_edge(n174, n173, tree.g).first;
  put(edge_weight_t(), tree.g, e, 58.61-38.64);
  
  Species* daubentonia = new Species("DaubentoniaMadagascariensis");
  vertex_descriptor vdaubentonia = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vdaubentonia, daubentonia);

  e = add_edge(n174, vdaubentonia, tree.g).first;
  put(edge_weight_t(), tree.g, e, 58.61);

  dmy = new Species("dmy144");
  vertex_descriptor n144 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n144, dmy);

  e = add_edge(n174, n144, tree.g).first;
  put(edge_weight_t(), tree.g, e, 68.74-58.61);
  
  Species* loris = new Species("LorisTardigradus");
  vertex_descriptor vloris = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vloris, loris);

  dmy = new Species("dmy177");
  vertex_descriptor n177 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n177, dmy);

  e = add_edge(n177, vloris, tree.g).first;
  put(edge_weight_t(), tree.g, e, 21.14);
  
  Species* nycticebus = new Species("NycticebusBengalensis");
  vertex_descriptor vnycticebus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vnycticebus, nycticebus);

  e = add_edge(n177, vnycticebus, tree.g).first;
  put(edge_weight_t(), tree.g, e, 21.14);

  dmy = new Species("dmy179");
  vertex_descriptor n179 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n179, dmy);

  e = add_edge(n177, n179, tree.g).first;
  put(edge_weight_t(), tree.g, e, 36.98-21.14);
  
  Species* perodicticus = new Species("PerodicticusPotto");
  vertex_descriptor vperodicticus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vperodicticus, perodicticus);

  e = add_edge(vperodicticus, n179, tree.g).first;
  put(edge_weight_t(), tree.g, e, 36.98);

  dmy = new Species("dmy184");
  vertex_descriptor n184 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n184, dmy);

  e = add_edge(n184, n179, tree.g).first;
  put(edge_weight_t(), tree.g, e, 40.34-36.98);
  
  Species* galago = new Species("GalagoSenegalensis");
  vertex_descriptor vgalago = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vgalago, galago);

  dmy = new Species("dmy182");
  vertex_descriptor n182 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n182, dmy);

  e = add_edge(n182, vgalago, tree.g).first;
  put(edge_weight_t(), tree.g, e, 15.37);
  
  Species* otolemur = new Species("OtolemurCrassicaudatus");
  vertex_descriptor votolemur = add_vertex(tree.g);
  put(Vertex_t(), tree.g, votolemur, otolemur);

  e = add_edge(n182, votolemur, tree.g).first;
  put(edge_weight_t(), tree.g, e, 15.37);

  e = add_edge(n184, n182, tree.g).first;
  put(edge_weight_t(), tree.g, e, 40.34-15.37);

  e = add_edge(n184, n144, tree.g).first;
  put(edge_weight_t(), tree.g, e, 68.74-40.34);

  e = add_edge(n185, n144, tree.g).first;
  put(edge_weight_t(), tree.g, e, 87.18-68.74);

  dmy = new Species("dmy187");
  vertex_descriptor n187 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n187, dmy);

  e = add_edge(n185, n187, tree.g).first;
  put(edge_weight_t(), tree.g, e, 92.26-87.18);
  
  Species* galeopterus = new Species("GaleopterusVariegatus(outgroup1)");
  vertex_descriptor vgaleopterus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vgaleopterus, galeopterus);

  e = add_edge(vgaleopterus, n187, tree.g).first;
  put(edge_weight_t(), tree.g, e, 92.26);

  dmy = new Species("dmy189");
  vertex_descriptor n189 = add_vertex(tree.g);
  put(Vertex_t(), tree.g, n189, dmy);

  e = add_edge(n189, n187, tree.g).first;
  put(edge_weight_t(), tree.g, e, 98.00-92.26);
  
  Species* tupaia = new Species("TupaiaBelangeri(outgroup2)");
  vertex_descriptor vtupaia = add_vertex(tree.g);
  put(Vertex_t(), tree.g, vtupaia, tupaia);

  e = add_edge(n189, vtupaia, tree.g).first;
  put(edge_weight_t(), tree.g, e, 98.00);
  
  Species* oryctolagus = new Species("OryctolagusCuniculus(outgroup3)");
  vertex_descriptor voryctolagus = add_vertex(tree.g);
  put(Vertex_t(), tree.g, voryctolagus, oryctolagus);

  e = add_edge(n189, voryctolagus, tree.g).first;
  put(edge_weight_t(), tree.g, e, 130.00);


  

  return tree;
}
