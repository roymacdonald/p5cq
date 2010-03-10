/*Constant Q efficient algorithm by Roy Macdonald macdonald.roy@gmail.com
based on the work by Judith C. Brown & Miller S. Puckette http://www.wellesley.edu/Physics/brown/pubs/effalgV92P2698-P2701.pdf
and the matlabimplementation by Benjamin Blankertz http://wwwmath1.uni-muenster.de/logik/org/staff/blankertz/constQ/constQ.html

*/
import krister.Ess.*;

AudioChannel myChannel;
ConstQ constQ;
float [] audio;
int sampRate;
void setup() {
  String loadPath = selectInput("Select audio file (.wav, .aiff, .mp3) ...");
  if (loadPath == null) {
    println("No file was selected...");
    stop();
  }
  else{
    size(1024,700); 
    smooth();

    Ess.start(this);
    myChannel=new AudioChannel(loadPath);
    audio = new float [myChannel.samples.length];
    arrayCopy(myChannel.samples, audio);
    sampRate = floor(myChannel.sampleRate);
    myChannel.destroy();
    constQ=new ConstQ(55, 7040, 24, sampRate, 0.0054 );//float minFreq, float maxFreq, int bins per octave, sampleRate, float threshold for sparseKernel(don't change unless you know what it does)
    background(0);
    
   image(constQ.processAndDraw(audio),0,0);
  }
}
void draw(){

}


public void stop() {
  Ess.stop();
  exit();
}

