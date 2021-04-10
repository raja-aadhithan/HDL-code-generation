// -------------------------------------------------------------
// 
// File Name: hdl_project_verilog\hdlsrc\project_detect\MATLAB_Function.v
// Created: 2021-02-13 17:18:24
// 
// Generated by MATLAB 9.9 and HDL Coder 3.17
// 
// -------------------------------------------------------------


// -------------------------------------------------------------
// 
// Module: MATLAB_Function
// Source Path: project_detect/Pulse Detector/Local Peak/MATLAB Function
// Hierarchy Level: 2
// 
// -------------------------------------------------------------

`timescale 1 ns / 1 ns

module MATLAB_Function
          (threshold,
           DataBuff_0,
           DataBuff_1,
           DataBuff_2,
           DataBuff_3,
           DataBuff_4,
           DataBuff_5,
           DataBuff_6,
           DataBuff_7,
           DataBuff_8,
           DataBuff_9,
           DataBuff_10,
           MidSample,
           detected);


  input   signed [17:0] threshold;  // sfix18_En11
  input   signed [17:0] DataBuff_0;  // sfix18_En11
  input   signed [17:0] DataBuff_1;  // sfix18_En11
  input   signed [17:0] DataBuff_2;  // sfix18_En11
  input   signed [17:0] DataBuff_3;  // sfix18_En11
  input   signed [17:0] DataBuff_4;  // sfix18_En11
  input   signed [17:0] DataBuff_5;  // sfix18_En11
  input   signed [17:0] DataBuff_6;  // sfix18_En11
  input   signed [17:0] DataBuff_7;  // sfix18_En11
  input   signed [17:0] DataBuff_8;  // sfix18_En11
  input   signed [17:0] DataBuff_9;  // sfix18_En11
  input   signed [17:0] DataBuff_10;  // sfix18_En11
  output  signed [17:0] MidSample;  // sfix18_En11
  output  detected;


  wire signed [17:0] DataBuff [0:10];  // sfix18_En11 [11]
  reg signed [17:0] MidSample_1;  // sfix18_En11
  reg  detected_1;
  reg signed [18:0] CompareOut [0:10];  // sfix19_En11 [11]
  reg  [0:10] c;  // boolean [11]
  reg  y;
  reg signed [18:0] extend_temp;  // sfix19_En11
  reg signed [31:0] k;  // int32
  reg signed [18:0] sub_cast [0:10];  // sfix19_En11 [11]


  assign DataBuff[0] = DataBuff_0;
  assign DataBuff[1] = DataBuff_1;
  assign DataBuff[2] = DataBuff_2;
  assign DataBuff[3] = DataBuff_3;
  assign DataBuff[4] = DataBuff_4;
  assign DataBuff[5] = DataBuff_5;
  assign DataBuff[6] = DataBuff_6;
  assign DataBuff[7] = DataBuff_7;
  assign DataBuff[8] = DataBuff_8;
  assign DataBuff[9] = DataBuff_9;
  assign DataBuff[10] = DataBuff_10;

  always @* begin
    // Hardware friendly implementation of peak finder
    //
    // Function inputs:
    // * WindowLen - non-tunable parameter defined under Simulink->Edit Data
    // * threhold  - input port (connected to constant)
    // * DataBuff  - input port (buffering done using Simulink block)
    //
    // Function outputs:
    // * "detected" is set when MidSample is local max
    // Compare each value in the window to the middle sample via subtraction
    extend_temp = {DataBuff[5][17], DataBuff[5]};
    // this is a vector
    // if all values in the result are negative and the middle sample is
    // greater than a threshold, it is a local max
    y = 1'b1;

    for(k = 32'sd0; k <= 32'sd10; k = k + 32'sd1) begin
      sub_cast[k] = {DataBuff[k][17], DataBuff[k]};
      CompareOut[k] = sub_cast[k] - extend_temp;
      c[k] = CompareOut[k] <= 19'sb0000000000000000000;
      y = y && c[k];
    end

    detected_1 = y && (DataBuff[5] > threshold);
    MidSample_1 = DataBuff[5];
  end



  assign MidSample = MidSample_1;

  assign detected = detected_1;

endmodule  // MATLAB_Function

