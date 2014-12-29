//Written by Nao Hiranuma
//Generates a raw sequence file using genbank files.
import java.util.*;
import java.io.*;

public class DataGenerator{
    ArrayList<Species> specieslist;
    Random rng;
    
    public DataGenerator(){
	rng = new Random();
        specieslist = new ArrayList<Species>();
        loadSpecies(loadfilenames());
    }

    public ArrayList<String> loadfilenames(){
        String path = "../GenbankData";
        String files;
        File folder = new File(path);
        File[] listOfFiles = folder.listFiles();
        
        ArrayList<String> filenames = new ArrayList<String>();
	int counter = 0;
        for(int i=0; i<listOfFiles.length; i++){
            if(listOfFiles[i].isFile()){
		counter++;
                filenames.add(listOfFiles[i].getName());
            }
        }
        System.out.println(counter+" species were found...");  
        return filenames;
    }
    
    public void loadSpecies(ArrayList<String> filenames){
        for(int i=0; i<filenames.size(); i++){
            specieslist.add(new Species(filenames.get(i)));
        }
    }
    
    public void writeresult(String writeTo, String geneName, int numS){
        try {
 
            File file = new File(writeTo);
            if (!file.exists()) {
                file.createNewFile();
            }
        
            FileWriter fw = new FileWriter(file.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fw);
	    double prob = 0;
	    double num_selected = 0;
	    int datasize = specieslist.size();
            for(int i=0; i<specieslist.size(); i++){
		prob = ((double)(numS - num_selected))/((double)(datasize - (i+1)));
		if(this.rng.nextFloat()<prob){
		    if(geneName.equals("cmp")){
			System.out.println("writing " + specieslist.get(i).name + " (" + specieslist.get(i).nucleotides.length()  +")");
			bw.write(specieslist.get(i).name + "\n");
			bw.write(specieslist.get(i).nucleotides + "\n");
		    }else{
			Gene g;
			for(int j=0; j<specieslist.get(i).genes.size();j++){
			    g = specieslist.get(i).genes.get(j);
			    if(g.genename.equals(geneName)){
				bw.write(specieslist.get(i).name + "\n");
				bw.write(specieslist.get(i).nucleotides.substring(g.genefrom, g.geneto)+ "\n");
				System.out.println("writing " + specieslist.get(i).name + "'s " + g.genename + " (" + (g.geneto-g.genefrom) +")");  
			    }
			}
		    }
		    num_selected ++;
                }
            }
        bw.close();
        }catch (IOException e) {
			e.printStackTrace();
		}
    }

    public void displayGenes(){
        Map<String, Integer> gcand = new HashMap<String, Integer>();
        ArrayList<Gene> gs;
        for(int i = 0; i < specieslist.size(); i++){
            for(int j = 0; j < specieslist.get(i).genes.size(); j++){
                gs = specieslist.get(i).genes;
                if(!gcand.containsKey(gs.get(j).genename)){
                    gcand.put(gs.get(j).genename, 1);
                }else{
                    gcand.put(gs.get(j).genename, gcand.get(gs.get(j).genename) + 1);
                }
            }
        }
        for(String n: gcand.keySet()){
            String key = n;
            Integer value = gcand.get(n);
            System.out.println( key + ": " + value + " species");
        }
    }
    
    public static void main(String [] args){
	if(args.length < 2){
	    System.out.println("Usage: java DataGenerator <std> <# of species> (<outfilename> <genename>)");
	    System.out.println("       for <std>");
	    System.out.println("            -y: takes gene name from stdout");
	    System.out.println("            -n: takes gene name from the command line arguments");
	    System.exit(1);
	}
	
	int numS = Integer.parseInt(args[0]);
	String flag = args[1];
	String outputfile;
	String geneName;

        DataGenerator dg = new DataGenerator();
	
	if(dg.specieslist.size() < numS){
	    System.out.println("The <# of species> you specified exceeds the datasize");
	    System.exit(1);
	}

        System.out.println("Species loaded\n");
        
	if(flag.equals("-y")){
	    System.out.println("Genes available:");
	    dg.displayGenes();
	    System.out.println();
	    Scanner ui = new Scanner(System.in);
	    System.out.print("Name your output file: ");
	    outputfile = ui.next();
	    System.out.print("Enter gene name (Say \"cmp\" for complete mtDNA sequences): ");
	    geneName = ui.next();
	}else{
	    outputfile = args[2];
	    geneName = args[3];
	}
        
        dg.writeresult(outputfile, geneName, numS);
        
        if(geneName.equals("cmp")){
        System.out.println("Complete mtDNA sequences were written to " + outputfile + ".");
        }else{
        System.out.println("DNA sequences for " + geneName + " were written to " + outputfile + ".");
        }
    }
}