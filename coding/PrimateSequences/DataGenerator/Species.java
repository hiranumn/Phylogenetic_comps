import java.util.*;
import java.io.*;

public class Species{
    String name;
    String nucleotides;
    ArrayList<Gene> genes;
    
    public Species(String filename){
        if(filename.contains(".")){
            name = filename.split("\\.")[0];
        }else{
            name = filename;
        }
        nucleotides = "";
        genes = new ArrayList<Gene>();
        loadData(filename);
	nucleotides = nucleotides.toUpperCase();
    }
    
    public void loadData(String filename){
        try {
	    File file = new File("../GenbankData/"+filename);
	    FileReader fileReader = new FileReader(file);
	    BufferedReader bufferedReader = new BufferedReader(fileReader);
	    String line;
            Boolean flag = false; //start loading
            String gene = "";
            int geneStart = 0;
            int geneEnd = 0;
            Boolean geneload = false;
			while ((line = bufferedReader.readLine()) != null) {
                if(flag == true){
                    line = line.replaceAll("\\s","");
                    line = line.replaceAll("/", "");
                    nucleotides += line.replaceAll("[0-9]","");
                }
                if(geneload == true){
                    Gene newGene = new Gene(line.split("\"")[1], geneStart, geneEnd);
                    genes.add(newGene);
                    geneload = false;
                }
                if(line.contains("ORIGIN")){
                    if(flag == false){
                        flag = true;
                    }
                }else if(line.length()>10){
                    if(line.substring(5,9).equals("gene")){
                        if(!line.contains("complement")){
                            line = line.split("\\s+")[2];
                            String[] info = line.split("\\.+");
                            geneStart = Integer.parseInt(info[0].replaceAll("[^\\d.]", ""));
                            geneEnd = Integer.parseInt(info[1].replaceAll("[^\\d.]", ""));
                            geneload = true;
                       }
                    }
                }

			}
			fileReader.close();
		}catch (IOException e) {
			e.printStackTrace();
		}
    }
}