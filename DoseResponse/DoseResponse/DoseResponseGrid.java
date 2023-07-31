package DoseResponse;

import HAL.Rand;
import HAL.Util;
import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GridWindow;
import static HAL.Util.*;

import java.util.ArrayList;
import java.util.List;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import HAL.Tools.FileIO;
import HAL.Gui.UIGrid;
import HAL.Interfaces.SerializableModel;

class Source extends AgentSQ2Dunstackable<DoseResponseGrid>{
    // Blood vessel source
    double conc;
    double xPos;
    double yPos;

    public Source(double conc, int x, int y) {
        super();
        Init(conc, x, y);
    }

    void Init(double conc,int x, int y){
        this.conc=conc;
        this.xPos=x;
        this.yPos=y;
    }

    void SetSource(PDEGrid2D diff){
        diff.Set(xPos,yPos,conc);
    }
}

class Cell extends AgentSQ2Dunstackable<DoseResponseGrid> {
    int genotype;
    double g_drugless;
    double ic50;
    double hillCoef;

    public int GetGenotype(){
        return this.genotype;
    }

    public double HillEqn(double conc){
        // convert conc to log10
        return this.g_drugless / (1 + Math.exp((this.ic50 - Math.log10(conc)) / this.hillCoef));
    }

    public double GetFitness(double conc){
        return HillEqn(conc);
    }
}
    

public class DoseResponseGrid extends AgentGrid2D<Cell> implements SerializableModel{
    
    // diffusion grid
    PDEGrid2D diff;
    Source src;
    double[] srcConc = new double[]{1000, 1000}; 
    int[] srcX = new int[]{50, 50};
    int[] srcY = new int[]{25, 75};
    // int[] srcX = new int[]{5,5};
    // int[] srcY = new int[]{3,7};
    Source[] srcList = new Source[srcX.length];
    double diffRate = 0.1;
    double diffBoundary = 0;
    double consumpRate = 0.5;
    
    // population initialization
    int xDim = 100;
    int yDim = 100; 
    int[] divHood = Util.VonNeumannHood(false);
    // public double initRadius = 10;
    public double initMutantProp = 0.1;
    String initGeometry = "circle"; // shape of initial population
    int initWidth = 100;
    double initDensity = 0.01;

    double dieProb = 0.1;
    double mutProb = 0.0001;
    double dt = 0.1; // time step in hours
    double[] GrowthRateList = new double[]{1.28949852, 1.14399848, 1.22802236, 0.93619847};
    double[] ic50List = new double[]{-0.49205992, 1.76224515,  1.39341393,  2.84653598};
    double hillCoef = -0.6824968;
    int n_genotype = 4;
    int nAllele = 2;
    int[] genotypeCounts = new int[n_genotype];

    // p.drugless_rates = [1.28949852, 1.14399848, 1.22802236, 0.93619847]
    // p.ic50 = [-0.49205992, 1.76224515,  1.39341393,  2.84653598]
    public int nReplicates = 1;
    int nTSteps = 1000;

    FileIO cellCountLogFile = null;
    String cellCountLogFileName = "./data/cellCountLog.csv";
    double logCellCountFrequency = 1;
    int tIdx = 0;
    boolean saveModelState = false;

    Rand rng = new Rand();
    int seed = 0;

    // Output - visualization
    UIGrid visCells;
    UIGrid visDiff;
    int scaleFactor = 2;
    boolean visualiseB = true;
    boolean saveFinalDiffImg = false;
    boolean saveFinalDiffGrid = false;
    boolean saveFinalPopGrid = false;
    String diffFileName = "~/tmp/diffGrid.csv";
    String popGridFileName = "~/tmp/popGrid.csv";
    public String imageOutDir = "./data/images/";
    int imageFrequency = -1;
    int pause = 1;
    int[] colorList = new int[]{RGB(0,1,0), RGB(0,0,1),
                                RGB(1,1,0), RGB(1,0,0)};

    // ------------------------------------------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------------------------------------------

    public DoseResponseGrid(int x, int y) {
        super(x, y, Cell.class);
        this.xDim = x;
        this.yDim = y;
    }

    public DoseResponseGrid(){
        super(100,100,Cell.class);
    }

    public DoseResponseGrid(int x, int y , double[] paramArr, double dt){
        super(x,y,Cell.class);
        // SetParameters(paramArr);
        this.dt = dt;
    }

    @Override
    public void SetupConstructors() {
        _PassAgentConstructor(Cell.class);
    }

    // ------------------------------------------------------------------------------------------------------------
    // Model setup and configuration
    // ------------------------------------------------------------------------------------------------------------

    public void SetSaveParams(String diffFileName, String cellCountLogFileName, 
                              String popGridFileName, boolean saveFinalDiffGrid,
                              boolean saveFinalPopGrid){
        this.diffFileName = diffFileName;
        this.cellCountLogFileName = cellCountLogFileName;
        this.saveFinalDiffGrid = saveFinalDiffGrid;
        this.saveFinalPopGrid = saveFinalPopGrid;
        this.popGridFileName = popGridFileName;
    }

    public void SetDiffParams(double[] srcConc, int[] srcX, int[] srcY, 
                              double diffRate, double diffBoundary, double consumpRate){
        this.srcConc = srcConc;
        this.srcX = srcX;
        this.srcY = srcY;
        this.diffRate = diffRate;
        this.diffBoundary = diffBoundary;
        this.srcList = new Source[srcX.length];
        this.consumpRate = consumpRate;
    }

    public void SetDiffParams(double[] srcConc,double diffRate, double consumpRate){
        this.srcConc = srcConc;
        this.diffRate = diffRate;
        this.srcList = new Source[srcX.length];
        this.consumpRate = consumpRate;
    }

    public void SetCellParams(double[] GrowthRateList, double[] ic50List, double hillCoef, 
                              int n_genotype, int nAllele, double mutProb, double dieProb){
        this.GrowthRateList = GrowthRateList;
        this.ic50List = ic50List;
        this.hillCoef = hillCoef;
        this.n_genotype = n_genotype;
        this.nAllele = nAllele;
        this.mutProb = mutProb;
        this.dieProb = dieProb;
    }

    public void SetCellParams(double mutProb, double dieProb){
        this.mutProb = mutProb;
        this.dieProb = dieProb;
    }
    /*
     * Set the initial population parameters. initGeometry is either "circle" or "square".
     */
    public void SetInitPopParams(String initGeometry, int initWidth, 
                                 double initDensity, double initMutantProp){
        this.initGeometry = initGeometry;
        // this.initRadius = initRadius;
        this.initWidth = initWidth;
        this.initDensity = initDensity;
        this.initMutantProp = initMutantProp;
    }

    public void SetSeed(int seed) {
        this.rng = new Rand(seed);
        // this.rn_ICs = new Rand(seed);
    }
    // public void ConfigureVisualisation(boolean visualiseB, int pause) {
    //     this.visualiseB = visualiseB;
    //     this.pause = pause;
    // }

    public void ConfigureImaging(boolean visualiseB, String imageOutDir, int imageFrequency,
                                 int scaleFactor, int pause, int[] colorList) {
        /*
        * Configure location of and frequency at which tumour is imaged.
        */
        this.visualiseB = visualiseB;
        this.imageOutDir = imageOutDir;
        this.imageFrequency = imageFrequency;
        this.scaleFactor = scaleFactor;
        this.pause = pause;
        this.colorList = colorList;
    }

    public void ConfigureImaging(boolean visualiseB, String imageOutDir, int imageFrequency,
                                 boolean saveFinalDiffImage) {
        /*
        * Configure location of and frequency at which tumour is imaged.
        */
        this.visualiseB = visualiseB;
        this.imageOutDir = imageOutDir;
        this.imageFrequency = imageFrequency;
        this.saveFinalDiffImg = saveFinalDiffImage;

    }

    public void ConfigureExperiment(int nReplicates, int nTSteps, double dt){
        this.nReplicates = nReplicates;
        this.nTSteps = nTSteps;
        this.dt = dt;
    }

    // ------------------------------------------------------------------------------------------------------------
    // Model initialization
    // ------------------------------------------------------------------------------------------------------------

    public void InitPopulation(){
        if ("circle".equals(this.initGeometry)) {
            InitPopCircle(this.initWidth, this.initMutantProp);
        }
        else if ("square".equals(this.initGeometry)) {
            InitPopSquare(this.initWidth, this.initMutantProp, this.initDensity);
        }
    }

    public void InitPopCircle(double radius, double mutantProb) {
        //get a list of indices that fill a circle at the center of the grid
        int[] tumorNeighborhood = CircleHood(true, radius);
        int hoodSize = MapHood(tumorNeighborhood, xDim / 2, yDim / 2);
        for (int i = 0; i < hoodSize; i++) {
            Cell c = NewAgentSQ(tumorNeighborhood[i]);
            if (rng.Double() < mutantProb) {
                int genotype = rng.Int(n_genotype-1) + 1; // any genotype other than 0
                c.genotype = genotype;
                c.g_drugless = GrowthRateList[genotype];
                c.ic50 = ic50List[genotype];
                c.hillCoef = hillCoef;
            } else {
                int genotype = 0;
                c.genotype = genotype;
                c.g_drugless = GrowthRateList[genotype];
                c.ic50 = ic50List[genotype];
                c.hillCoef = hillCoef;
            }
            genotypeCounts[c.genotype] += 1;
        }
    }

    public void InitPopSquare(int x, double mutantProb, double density){
        
        int hoodSize = x*x;
        for (int i = 0; i < hoodSize; i++){
            if (rng.Double() < density){
                Cell c = NewAgentSQ(i);
                if (rng.Double() < mutantProb) {
                    int genotype = rng.Int(n_genotype-1) + 1; // any genotype other than 0
                    c.genotype = genotype;
                    c.g_drugless = GrowthRateList[genotype];
                    c.ic50 = ic50List[genotype];
                    c.hillCoef = hillCoef;
                } else {
                    int genotype = 0;
                    c.genotype = genotype;
                    c.g_drugless = GrowthRateList[genotype];
                    c.ic50 = ic50List[genotype];
                    c.hillCoef = hillCoef;
                }
                genotypeCounts[c.genotype] += 1;
            }
            
        }
    }

    public void InitPDEGrid() {
        diff = new PDEGrid2D(xDim, yDim);
        // initialize sources
        for (int i = 0; i < srcX.length; i++) {
            Source src = new Source(srcConc[i], srcX[i], srcY[i]);
            srcList[i] = src;
        }
    }

    // ------------------------------------------------------------------------------------------------------------
    // Model step
    // ------------------------------------------------------------------------------------------------------------

    public void StepPDEGrid() {
        for (Source src : srcList) {
            src.SetSource(diff); // sets the concentration at the source's x,y position
        }
        diff.Diffusion(diffRate,diffBoundary);
        // field consumption
        // for (int i = 0; i < diff.length; i++) {
        //     diff.Mul(i, -consumpRate);
        // }
        diff.MulAll(-consumpRate);
        diff.Update();
    }

    public List<Integer> genNeighbors(int genotype, int nAllele) {
        List<Integer> neighbors = new ArrayList<>();
        int mask = (1 << nAllele) - 1;

        for (int m = 0; m < nAllele; m++) {
            int neighbor = genotype ^ (1 << m);
            neighbors.add(neighbor & mask);
        }

        return neighbors;
    }

    public void StepModel(){
        StepCells();
        StepPDEGrid();
    }

    public void StepCells() {
        int currPos;

        for (Cell cell : this) {
            currPos = cell.Isq();
            // cell.StepCell(dieProb);
            if (rng.Double() < dieProb * dt){
                cell.Dispose();
                genotypeCounts[cell.genotype] -= 1;
            }
            else{
                // get the drug concentration
                double currConc = diff.Get(currPos);
                double divProb = cell.GetFitness(currConc) * dt;

                // System.out.print(currConc + " ");
                // System.out.println(divProb);

                if (rng.Double() < divProb) {
                    int options = MapEmptyHood(divHood, cell.Xsq(), cell.Ysq());
                    if (options > 0) {

                        int iDaughter = divHood[rng.Int(options)];
                        Cell daughter = NewAgentSQ(iDaughter);
                        
                        if (rng.Double() < mutProb){
                            // mutate
                            List<Integer> neighbors = genNeighbors(cell.genotype, nAllele);
                            int newGenotype = neighbors.get(rng.Int(neighbors.size()));
                            daughter.genotype = newGenotype;
                            daughter.g_drugless = GrowthRateList[newGenotype];
                            daughter.ic50 = ic50List[newGenotype];
                            daughter.hillCoef = hillCoef;
                            
                        }
                        else{
                            // inherit properties from parent
                            daughter.genotype = cell.genotype;
                            daughter.g_drugless = cell.g_drugless;
                            daughter.ic50 = cell.ic50;
                            daughter.hillCoef = cell.hillCoef;
                            
                        }
                        genotypeCounts[daughter.genotype] += 1;
                    }
                }
            }
        }
    }
    
    public static void main(String[] args) {

        int x=100;
        int y=100;
        int timesteps=1000;
        GridWindow visCells=new GridWindow(x,y,5,true,null,true);
        GridWindow visDiff=new GridWindow(x,y,5,true,null,true);
        // double dieProb=0.01;
    
        // GridWindow win = new GridWindow(x,y,10);
        DoseResponseGrid model = new DoseResponseGrid(x,y);

        model.InitPopulation();
        model.InitPDEGrid();

        // get the number of resistant cells in the tumor
    
        for (int i = 0; i < timesteps; i++) {
            // model step
            model.StepModel();
            model.Draw(visCells,visDiff);
            visDiff.TickPause(1);
        }
        // model.SaveDiffGrid(model.diff, "/Users/eshanking/repos/HAL_fds/DoseResponse/data/diffGrid.csv");
    }

    public void Run() throws IOException {

        UIGrid currVis = new UIGrid(xDim, yDim, scaleFactor, visualiseB); // For head-less run
        this.visCells = currVis;

        Boolean completedSimulationB = false;
        Boolean logged = false;
        // System.out.println(this.initGeometry);
        InitPopulation();
        // initialize the diffusion grid
        InitPDEGrid();
        if (cellCountLogFile==null && cellCountLogFileName!=null) {InitialiseCellLog(this.cellCountLogFileName);}

        SaveCurrentCellCount(0);
        SaveTumourImage(tIdx);
        
        completedSimulationB = false;
        while (!completedSimulationB) {
            // currVis.TickPause(pause);
            Draw(currVis);
            StepModel();
            logged = SaveCurrentCellCount(tIdx);
            SaveTumourImage(tIdx);
            tIdx++;
            completedSimulationB = (tIdx>nTSteps)?true:false;
        }
        if (saveFinalDiffImg){
            SaveDiffImage(tIdx);
        }
        if (saveFinalDiffGrid){
            SaveDiffGrid(diffFileName);
        }

        if (saveFinalPopGrid){
            SavePopGrid(popGridFileName);
        }
        this.Close(logged);
    }
    // ------------------------------------------------------------------------------------------------------------
    // Manage and save output
    // ------------------------------------------------------------------------------------------------------------

    public void InitialiseCellLog(String cellCountLogFileName) {
        cellCountLogFile = new FileIO(cellCountLogFileName, "w");
        WriteLogFileHeader();
        this.cellCountLogFileName = cellCountLogFileName;
        this.logCellCountFrequency = 1;
    }

    private void WriteLogFileHeader() {
        // cellCountLogFile.Write("TIdx,Time,NCells_S,NCells_R,NCells,DrugConcentration,rS,rR,mS,mR,dS,dR,dD_div_S,dD_div_R,dt");

        // cellCountLogFile.Write("TIdx,Time,NCells_S,NCells_R,NCells");
        cellCountLogFile.Write("TIdx,Time,Pop,n_gen0,n_gen1,n_gen2,n_gen3,mutProb,dieProb,diffRate,consumpRate,src1_Conc,src2_Conc,src1_X,src2_X,src1_Y,src2_Y");
        
        cellCountLogFile.Write("\n");
    }
    public void Close() {
        if (cellCountLogFile!=null) {cellCountLogFile.Close();}
    }

    public Boolean SaveCurrentCellCount(int currTimeIdx) {
        Boolean successfulLog = false;
        if ((currTimeIdx % (int) (logCellCountFrequency)) == 0 && logCellCountFrequency > 0) {
            cellCountLogFile.WriteDelimit(GetModelState(),",");
            // if (extraSimulationInfoNames!=null) {
            //     cellCountLogFile.Write(",");
            //     cellCountLogFile.WriteDelimit(extraSimulationInfo, ",");
            // }
            cellCountLogFile.Write("\n");
            successfulLog = true;
        }
        return successfulLog;
    }

    public double[] GetModelState() {
   
        return new double[] {tIdx, tIdx*dt, Pop(), genotypeCounts[0], genotypeCounts[1], 
                             genotypeCounts[2], genotypeCounts[3], mutProb, dieProb, diffRate,
                             consumpRate, srcConc[0], srcConc[1], srcX[0], srcX[1], srcY[0], srcY[1]};
        
    }

    public void SaveTumourImage(int currTimeIdx) {
        if (imageFrequency > 0 && (currTimeIdx % (int) (imageFrequency/dt)) == 0) {
            this.visCells.ToPNG(imageOutDir +"/img_t_"+currTimeIdx*dt+".png");
        }
    }

    public void SaveDiffImage(int tIdx) {

        if (this.visDiff == null){
            this.visDiff = new UIGrid(xDim, yDim, scaleFactor, visualiseB);
        }
        Draw(this.visCells,this.visDiff);

        this.visDiff.ToPNG(imageOutDir +"/diffGrid_"+tIdx*dt+".png");

    }

    public void SaveDiffGrid(String diffFileName) throws IOException {
        try{
            File diffGridFile = new File(diffFileName);
            diffGridFile.createNewFile();
            FileWriter diffGridFileWriter = new FileWriter(diffFileName);
            double[] field = diff.GetField();
            
            diffGridFileWriter.write("x,y,field\n");
            // int lineNum = 0;
            for (int i = 0; i < field.length; i++) {
                diffGridFileWriter.write(i % xDim + "," + i / xDim + "," + field[i] + "\n");
                // lineNum++;
            }
            diffGridFileWriter.close();
            // System.out.println("Wrote "+lineNum+" lines to "+diffFileName);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /*
     * Save the population grid to a file
     */
    public void SavePopGrid(String popFileName) throws IOException {
        try{
            File popGridFile = new File(popFileName);
            popGridFile.createNewFile();
            FileWriter popGridFileWriter = new FileWriter(popFileName);
            
            popGridFileWriter.write("x,y,genotype\n");
            int currPos;
            for (Cell cell : this){
                currPos = cell.Isq();
                popGridFileWriter.write(currPos % xDim + "," + currPos / xDim + "," + cell.GetGenotype() + "\n");
            }
            popGridFileWriter.close();
            // System.out.println("Wrote "+lineNum+" lines to "+diffFileName);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    // public void SaveDiffGrid(String diffFileName) {
    //     // get the diffusion grid values
    //     double[] field = diff.GetField();
    //     // new fileio object
    //     FileIO diffGridFile = new FileIO(diffFileName, "w");
    //     // write the header
    //     diffGridFile.Write("x,y,field\n");
    //     // write the values
    //     int lineNum = 0;
    //     for (int i = 0; i < field.length; i++) {
    //         System.out.println(i % xDim + "," + i / xDim + "," + field[i]);
    //         diffGridFile.Write(i % xDim + "," + i / xDim + "," + field[i] + "\n");
    //         lineNum++;
    //     }
    //     System.out.println("Wrote "+lineNum+" lines to "+diffFileName);

    // }
    public void SaveDiffGrid(PDEGrid2D diffGrid, String diffFileName) {
        double[] field = diffGrid.GetField();
        FileIO diffGridFile = new FileIO(diffFileName, "w");
        diffGridFile.Write("x,y,field\n");
        
        for (int i = 0; i < field.length; i++) {
            diffGridFile.Write(i % xDim + "," + i / xDim + "," + field[i] + "\n");
        }
    }

    public void Close(Boolean logged) {
        if (!logged) {
            tIdx--;
            SaveCurrentCellCount(0);}
        if (cellCountLogFile!=null) {cellCountLogFile.Close();}
    }
    public void SaveModelState(String stateFileName) {
        // Can't have active pointers when saving the model. So, close everything here.
        if (cellCountLogFile!=null) {
            cellCountLogFile.Close();
            cellCountLogFile = null;
        }
        if (visCells!=null) {visCells = null;}
        SaveState(this,stateFileName);
    }
    /*
     * Draw cells and diffusion grid
     */
    public void Draw(UIGrid visCells, UIGrid visDiff){
        for (Cell cell : this) {
            visCells.SetPix(cell.Isq(),colorList[cell.genotype]);//draw sources and sinks
        }
        for (int i = 0; i < length; i++) {//length of the Grid
            //visDiff.SetPix(i,SetAlpha(HeatMapRGB(diff.Get(i)*4),diff.Get(i)*4));
            visDiff.SetPix(i,HeatMapRGB(diff.Get(i)*4));
        }
    }
    /*
     * Draw cells only
     */
    public void Draw(UIGrid visCells){
        for (Cell cell : this) {
            visCells.SetPix(cell.Isq(),colorList[cell.genotype]);//draw sources and sinks
        }
    }
}
