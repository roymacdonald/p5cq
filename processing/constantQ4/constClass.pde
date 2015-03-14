class ConstQ {
  SparseArray [] sparKernel;
  float Q;
  int K;
  int fftLen;
  int sparLen;
  //////////////////////////////////////----Const Q CONSTRUCTOR---/////////////////////////////////////////////////////////////////////
  ///
  ConstQ(float minFreq, float maxFreq, int bins, int fs, float thresh){
    sparseKernel(minFreq, maxFreq, bins, fs, thresh);
  }
  //////////////////////////////////////----Sparse Kernel----/////////////////////////////////////////////////////////////////////////
  private void sparseKernel(float minFreq, float maxFreq, int bins, int fs, float thresh){

    Q= 1/(pow(2,1.0f/bins)-1);                                         //Q= 1/(2^(1/bins)-1);

    K= ceil(bins*log2(maxFreq/minFreq));                                 //K= ceil( bins * log2(maxFreq/minFreq) );
    println("Q: "+Q);
    println("K: "+K);

    fftLen= round(pow(2, nextpow2(ceil((Q*fs)/minFreq))));               //fftLen= 2^nextpow2( ceil(Q*fs/minFreq) );

    println("fftLen: "+fftLen);


    float [][] tempKernel= new float [fftLen][2];

    Complex [] specKernel = new Complex[fftLen];       //tempKernel= zeros(fftLen, 1);
    sparKernel=new SparseArray[K];
    for(int i=0; i<fftLen;i++){//////////////  Fill tempKernel with zeros.   
      tempKernel[i][0]=0;//real
      tempKernel[i][1]=0;//img
    }
    //---------------------
    FFTX fft1= new FFTX(fftLen);

    float yy=0;


    for(int k=K;k>=1;k--){            //for k= K:-1:1;
      int len= ceil(Q*fs/(minFreq*pow(2,((k-1)/(float)bins))));                     //len= ceil( Q * fs / (minFreq*2^((k-1)/bins)) );
      float [] hamm = hamming(len);                                          // tempKernel(1:len)= hamming(len)/len .* exp(2*pi*i*Q*(0:len-1)'/len);
      float [] tmp= new float [len];
      float u = 2*PI*Q/len;
      for(int n=0;n<len;n++){
        float hl =hamm[n]/len;
        float un=u*n;
        tempKernel[n][0] = hl * cos(un);//real
        tempKernel[n][1] = hl * sin(un);//img
        tmp[n] = tempKernel[n][0];
      }
      //---------------------
      fft1.forward(tempKernel);                                              //specKernel= fft(tempKernel);

      for(int c=0; c<specKernel.length; c++){
        specKernel[c]=new Complex(fft1.real[c],fft1.imag[c]);
      }                                                                      //---------------------
      sparKernel[k-1]= new SparseArray(); 
      for(int z=0; z<specKernel.length; z++){                                //specKernel(find(abs(specKernel)<=thresh))= 0;
        if(specKernel[z].magnitude()>=thresh){
          sparKernel[k-1].addSparse(new Complex(specKernel[z].r/=fftLen,specKernel[z].i*=-fftLen),z);

        }
      }                                                                     //---------------------


      println(k);
    }


    //---------------------
    sparLen =  round(pow(2,nextpow2(sparseSize())))*2;
    println("sparKernel length: "+sparLen);
    
  }

  //////////////////////////////////////----Guardar----////////////////////////////////////////////////////////////////////////////////
  void guardar(){
    String [] st = new String [K]; 
    for(int c=0; c<K; c++){
      String ss = "";  
      for(int i=0; i<sparKernel[c].sparse.length; i++){
        ss+= str(sparKernel[c].sparse[i].complex.magnitude());
        ss+= "  ";
      }
      st[c]=ss;
    }
    saveStrings("sparKernel.txt", st);

  }
  //////////////////////////////////////----Dibujar Sparse----/////////////////////////////////////////////////////////////////////////
  /*Just for debugging. It draws the sparse kernel like the one showed in the paper by Brown && Puckette
  
  */
  PImage dibujarSparse(){
    PGraphics pg = createGraphics(width, height,P3D);
    float factor = 6*width/float(fftLen);
    println("factor x: "+factor);
    pg.beginDraw();
    pg.background(0);
    pg.stroke(255);
    pg.smooth();
    float maxReal;
    float yPos=0;
    for(int c=0; c<K; c++){
      float [] real=new float [sparKernel[c].sparse.length];
      for(int i= 0; i<sparKernel[c].sparse.length; i++){
        real[i]=sparKernel[c].sparse[i].r;
      }
      maxReal =(height/float(K*2))/max(real);
      if(c>0){
        yPos=c*height/float(K);
      }
      println("yPos: "+yPos+"  max: "+maxReal);
      for(int i=1; i<sparKernel[c].sparse.length; i++){
        float x1 = factor*sparKernel[c].sparse[i-1].index;
        float y1 = (sparKernel[c].sparse[i-1].r*maxReal)+yPos;
        float x2 = factor*sparKernel[c].sparse[i].index;
        float y2 = (sparKernel[c].sparse[i].r*maxReal)+yPos;
        pg.line(x1,y1,x2,y2);
      }

    }
    pg.endDraw();
    return pg;  
  }
  int sparseSize(){
    int m=0;
    for(int i = 0; i<sparKernel.length;i++){
      for(int j=0; j<sparKernel[i].sparse.length;j++){
        if(sparKernel[i].sparse[j].index>m){
          m=sparKernel[i].sparse[j].index;
        }
      }
    }
    return m;
  }
  //////////////////////////////////////----Process & draw----///////////////////////////////////////////////////////////////////////
  /*
  So far this only draws the complete audio file processed by the cosntant q transform.
  The complete audiofile is chopped into the horizontal number of pixels, this might yield into having chunks that share sample values. This is helpful when viewing shot audio files, so to have an "increased resolution"
  Feel free to modify this part.
  */
  PImage processAndDraw(float[] a){
    PGraphics pg = createGraphics(width, height);//,P3D);
    float factorH= height/float(K);
    //int cuantos = floor(a.length/float(fftLen));
    int cuantos =width-1;
    float factorW = width/float(cuantos);
    int inc = floor(a.length/float(width));
    //int inc = 100;
    //int inc = sparLen;
    //int inc =400;
    println("sample Increment: "+inc);
    pg.beginDraw();
    pg.background(0);
    pg.noStroke();
    // pg.smooth();
    for(int i=0; i<cuantos && i<width; i++){
      float [] temp;
      if(a.length-(i*inc)>=sparLen){
        temp = process(subset(a,i*inc,sparLen));
      }
      else{
        temp = process(subset(a,i*inc,a.length-i*inc));
      }
      for(int j =1; j<K;j++){
        pg.stroke(temp[j]/float(10000));
        pg.line(i, factorH*(j-1), i, factorH*(j));
      }
      //println("Procesando Bin "+ i+ ":" +cuantos);
    }
    pg.endDraw();
    return pg;  
  }




//////////////////////////////////////----process----//////////////////////////////////////////////////////////////////////////

  float [] process(float [] x){// x must be a row vector  
    FFTX fft = new FFTX(sparLen);
    if(x.length>sparLen){
      x=subset(x, 0, sparLen);
    }
    else if(x.length<sparLen){
      int largozeros=sparLen-x.length;
      float[]zeros = new float[largozeros];
      for(int z=0; z<largozeros;z++){
        zeros[z]=0;
      }
      x=concat(x,zeros);
    }
    fft.forward(x);
    float [] temp = new float[sparKernel.length];
    for(int i=0; i<sparKernel.length; i++){
      for(int j=0; j<sparKernel[i].sparse.length; j++){
        temp[i]+=fft.getSpectrum()[sparKernel[i].sparse[j].index]*sparKernel[i].sparse[j].complex.magnitude();
      }
    }
    //  println(temp);
    return temp;
  }


  //////////////////////////////////////----Next power of 2----////////////////////////////////////////////////////////////////////////

  int nextpow2(int in){
    int i=1;
    while(true){
      i++;

      if(pow(2,i)>=in){
        return i;
      }
    }
  }
//////////////////////////////////////----Complex array to SparseComplex array function----//////////////////////////////////////////
Complex [] sparse(Complex [] c){
  Complex [] spar={
  };
  for(int a=0; a<c.length; a++){
    if(c[a].r !=0 && c[a].i!=0){
      spar =(Complex[]) append(spar, c[a]);  
    }
  }
  return spar;
}
//////////////////////////////////////----Log base 2----/////////////////////////////////////////////////////////////////////////////
  float log2 (float x) {
    return (log(x) / log(2));
  }
//////////////////////////////////////----Hamming window of length n----/////////////////////////////////////////////////////////////
  float [] hamming(int n){
    float [] temp=new float [n];
    for(int i = 0; i<n;i++){
      temp[i]= 0.53836-(0.46164*cos(TWO_PI*i/(n-1)));
    }
    return temp;
  }
}

//////////////////////////////////////----SparseArray Class----//////////////////////////////////////////////////////////////////////
  class SparseArray{
    SparseComplex [] sparse={
    };  
    SparseArray(){
    }
    void addSparse(Complex c, int ind){
      sparse=(SparseComplex[])append(sparse, new SparseComplex(c,ind));   
    }
  }
//////////////////////////////////////----SparseComplex Class----////////////////////////////////////////////////////////////////////
  class SparseComplex{
    int index;
    float r;
    float i;
    Complex complex;
    SparseComplex (Complex c,int ind){
      complex =c;
      index=ind;
      r=c.r;
      i=c.i;
    }
  } 

//////////////////////////////////////----Complex number class----/////////////////////////////////////////////////////////////
class Complex{
  float r, i;
  Complex(float _r, float _i){
    r=_r;
    i=_i;
  }
  void zero(){
    r=0;
    i=0;
  }
  float magnitude(){
    return (float)Math.sqrt(r*r + i*i);
  }
}

Complex mlt(Complex a, Complex c){ //(a + bi)(c + di) = (ac − bd) + (bc + ad)i
  Complex res=new Complex(a.r*c.r - a.i*c.i, a.i*c.r + a.r*c.i);
  return res;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//         FFT
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Standard Fast Fourier Transform (fft) implementation
*/
class FFTX{
  FFTX(int ts) {
    timeSize = ts;
    // sampleRate = (int)sr;
    //bandWidth = (2.0 / (float)timeSize) * ((float)sampleRate / 2.0);
    allocateArrays();
    if((timeSize & timeSize - 1) != 0) {
      throw new IllegalArgumentException("FFT: timeSize must be a power of two.");
    }
    else {
      rvrs = new int[timeSize];
      sinlookup = new float[timeSize];
      coslookup = new float[timeSize];
      buildReverseTable();
      buildTrigTables();
      return;
    }
  }

  protected void allocateArrays() {
    spectrum = new float[timeSize / 2 + 1];
    real = new float[timeSize];
    imag = new float[timeSize];
  }

  private void fft() {
    for(int halfSize = 1; halfSize < real.length; halfSize *= 2) {
      float phaseShiftStepR = coS(halfSize);
      float phaseShiftStepI = siN(halfSize);
      float currentPhaseShiftR = 1.0;
      float currentPhaseShiftI = 0.0;
      for(int fftStep = 0; fftStep < halfSize; fftStep++) {
        for(int i = fftStep; i < real.length; i += 2 * halfSize) {
          int off = i + halfSize;
          float tr = currentPhaseShiftR * real[off] - currentPhaseShiftI * imag[off];
          float ti = currentPhaseShiftR * imag[off] + currentPhaseShiftI * real[off];
          real[off] = real[i] - tr;
          imag[off] = imag[i] - ti;
          real[i] += tr;
          imag[i] += ti;
        }
        float tmpR = currentPhaseShiftR;
        currentPhaseShiftR = tmpR * phaseShiftStepR - currentPhaseShiftI * phaseShiftStepI;
        currentPhaseShiftI = tmpR * phaseShiftStepI + currentPhaseShiftI * phaseShiftStepR;
      }
    }
  }

  public void forward(float [][] buffer) {

    bitReverseSamples(buffer);
    fft();
    fillSpectrum();
    return;  
  }
  public void forward(float [] buffer) {
    if(buffer.length != timeSize) {
      println("FFT.forward: The length of the passed sample buffer must be equal to timeSize().");
      return;
    } 
    else {
      bitReverseSamples(buffer);
      fft();
      fillSpectrum();
      return;
    }
  }
  private void buildReverseTable() {
    int N = timeSize;
    rvrs[0] = 0;
    int lmt = 1;
    //println("N: "+N+"  N/2: "+ N/ 2);
    for(int bit = N/ 2; lmt < N; bit >>= 1) {
      for(int i = 0; i < lmt; i++)
        rvrs[i + lmt] = rvrs[i] + bit;
      lmt <<= 1;
    }
  }

  private void bitReverseSamples(float [][] samples) {
    for(int i = 0; i < timeSize; i++) {
      real[i] = samples[rvrs[i]][0];
      imag[i] = samples[rvrs[i]][1];
    }
  }
  private void bitReverseSamples(float [] samples) {
    for(int i = 0; i < timeSize; i++) {
      real[i] = samples[rvrs[i]];
    }
  }
  private float siN(int i) {
    return sinlookup[i];
  }

  private float coS(int i) {
    return coslookup[i];
  }

  private void buildTrigTables() {
    int N = timeSize;
    for(int i = 0; i < N; i++) {
      sinlookup[i] = sin(-PI /float(i));
      coslookup[i] = cos(-PI /float(i));
    }
  }
  protected void fillSpectrum() {
    for(int i = 0; i < spectrum.length; i++)
      spectrum[i] = (float)Math.sqrt(real[i] * real[i] + imag[i] * imag[i]);
  }
  public float [] getSpectrum(){
    return spectrum;
  }

  private int rvrs[];
  private float sinlookup[];
  private float coslookup[];
  protected int timeSize;
  protected float real[];
  protected float imag[];
  protected float spectrum[];
}







