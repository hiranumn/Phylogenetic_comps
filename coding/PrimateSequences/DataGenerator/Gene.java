public class Gene{
    String genename;
    int genefrom;
    int geneto;
    
    public Gene(String gn, int gf, int gt){
        genename = gn;
        genefrom = gf;
        geneto= gt;
    }

    public String toString(){
        return genename+": "+genefrom+" .. "+ geneto;// + " size: " + 
    }
}
