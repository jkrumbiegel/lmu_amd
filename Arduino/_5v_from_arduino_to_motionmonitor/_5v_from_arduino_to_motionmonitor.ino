// 05.07.2016
// Sketch to send a 5v pulse to the MC USB-1616HS-BNC device to synchronize the MotionMonitor
// "p" sent from Matlab sends one pulse of 50ms

void setup() {
  Serial.begin(9600);
  pinMode(13,OUTPUT); // 5v will be applied to this pin on HIGH state
  digitalWrite(13,LOW); //make sure pin is at 0v
}

void loop() {
  if(Serial.available() > 0) {
    int inByte = Serial.read();
    if (inByte == 'p') {
      digitalWrite(13,HIGH);  // set pin 13 to 5v
      delay(50);              // wait 50 ms
      digitalWrite(13,LOW);   // reset pin 13 to 0v
    }    
  }
}
