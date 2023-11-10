/*--------------------------------------------------------------------------------------------------------
This library aims to provide Sentinel-2 biophysical parameter retrievals through GEE, based on the 
S2ToolBox methodology.
For algorithm details, see the original ATBD: https://step.esa.int/docs/extra/ATBD_S2ToolBox_L2B_V1.1.pdf

Currently, only FAPAR and LAI have been ported to GEE. fCOVER can be done as well.
Input should always be Sentinel-2 L2A products. 

This is a lot of neural network parameters, and there has been --no-- thorough validation of this code.
Please use at your own risk and provide feedback to:

kristofvantricht@gmail.com

--------------------------------------------------------------------------------------------------------
*/
var degToRad = ee.Image(Math.PI / 180);


//--------------------
// FAPAR 8-BAND
//--------------------

exports.get_fapar = function(img) {
  img = descale(img, 0.0001);
  var inputoormask = maskinputOOR(img);
  var b03_norm = normalize(img.select('B3'), 0, 0.253061520471542);
  var b04_norm = normalize(img.select('B4'), 0, 0.290393577911328);
  var b05_norm = normalize(img.select('B5'), 0, 0.305398915248555);
  var b06_norm = normalize(img.select('B6'), 0.006637972542253, 0.608900395797889);
  var b07_norm = normalize(img.select('B7'), 0.013972727018939, 0.753827384322927);
  var b8a_norm = normalize(img.select('B8A'), 0.026690138082061, 0.782011770669178);
  var b11_norm = normalize(img.select('B11'), 0.016388074192258, 0.493761397883092);
  var b12_norm = normalize(img.select('B12'), 0, 0.493025984460231);
  var viewZen_norm = normalize(ee.Image(img.getNumber('MEAN_INCIDENCE_ZENITH_ANGLE_B8')).multiply(degToRad).cos(), 0.918595400582046, 1);
  var sunZen_norm  = normalize(ee.Image(img.getNumber('MEAN_SOLAR_ZENITH_ANGLE')).multiply(degToRad).cos(), 0.342022871159208, 0.936206429175402);
  var relAzim_norm = (ee.Image(img.getNumber('MEAN_SOLAR_AZIMUTH_ANGLE')).subtract(ee.Image(img.getNumber('MEAN_INCIDENCE_AZIMUTH_ANGLE_B8'))).multiply(degToRad)).cos()

  var n1 = neuron1_fapar8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm);
  var n2 = neuron2_fapar8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm);
  var n3 = neuron3_fapar8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm);
  var n4 = neuron4_fapar8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm);
  var n5 = neuron5_fapar8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm);
  
  var l2 = layer2_fapar8(n1, n2, n3, n4, n5);
  
  var fapar = denormalize(l2, 0.000153013463222, 0.977135096979553);
  return fapar.updateMask(inputoormask)
}


function neuron1_fapar8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm) {
  var sum =
	ee.Image(-0.887068364040280)
	.add(ee.Image(0.268714454733421).multiply(b03_norm))
	.subtract(ee.Image(0.205473108029835).multiply(b04_norm))
	.add(ee.Image(0.281765694196018).multiply(b05_norm))
	.add(ee.Image(1.337443412255980).multiply(b06_norm))
	.add(ee.Image(0.390319212938497).multiply(b07_norm))
	.subtract(ee.Image(3.612714342203350 ).multiply(b8a_norm))
	.add(ee.Image(0.222530960987244).multiply(b11_norm))
	.add(ee.Image(0.821790549667255).multiply(b12_norm))
	.subtract(ee.Image(0.093664567310731).multiply(viewZen_norm))
	.add(ee.Image(0.019290146147447).multiply(sunZen_norm))
	.add(ee.Image(0.037364446377188).multiply(relAzim_norm));

  return tansig(sum);
  
}

function neuron2_fapar8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm) {
  var sum =
	ee.Image(0.320126471197199)
	.subtract(ee.Image(0.248998054599707).multiply(b03_norm))
	.subtract(ee.Image(0.571461305473124).multiply(b04_norm))
	.subtract(ee.Image(0.369957603466673).multiply(b05_norm))
	.add(ee.Image(0.246031694650909).multiply(b06_norm))
	.add(ee.Image(0.332536215252841).multiply(b07_norm))
	.add(ee.Image(0.438269896208887).multiply(b8a_norm))
	.add(ee.Image(0.819000551890450).multiply(b11_norm))
	.subtract(ee.Image(0.934931499059310).multiply(b12_norm))
	.add(ee.Image(0.082716247651866).multiply(viewZen_norm))
	.subtract(ee.Image(0.286978634108328).multiply(sunZen_norm))
	.subtract(ee.Image(0.035890968351662).multiply(relAzim_norm));
	
  return tansig(sum);
}

function neuron3_fapar8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm) {
  var sum =
	ee.Image(0.610523702500117)
	.subtract(ee.Image(0.164063575315880).multiply(b03_norm))
	.subtract(ee.Image(0.126303285737763).multiply(b04_norm))
	.subtract(ee.Image(0.253670784366822).multiply(b05_norm))
	.subtract(ee.Image(0.321162835049381).multiply(b06_norm))
	.add(ee.Image(0.067082287973580).multiply(b07_norm))
	.add(ee.Image(2.029832288655260).multiply(b8a_norm))
	.subtract(ee.Image(0.023141228827722).multiply(b11_norm))
	.subtract(ee.Image(0.553176625657559).multiply(b12_norm))
	.add(ee.Image(0.059285451897783).multiply(viewZen_norm))
  .subtract(ee.Image(0.034334454541432).multiply(sunZen_norm))
	.subtract(ee.Image(0.031776704097009).multiply(relAzim_norm));
	
  return tansig(sum);
}

function neuron4_fapar8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm) {
  var sum =
	ee.Image(-0.379156190833946)
	.add(ee.Image(0.130240753003835).multiply(b03_norm))
	.add(ee.Image(0.236781035723321).multiply(b04_norm))
	.add(ee.Image(0.131811664093253).multiply(b05_norm))
	.subtract(ee.Image(0.250181799267664).multiply(b06_norm))
	.subtract(ee.Image(0.011364149953286).multiply(b07_norm))
	.subtract(ee.Image(1.857573214633520).multiply(b8a_norm))
	.subtract(ee.Image(0.146860751013916).multiply(b11_norm))
	.add(ee.Image(0.528008831372352).multiply(b12_norm))
	.subtract(ee.Image(0.046230769098303).multiply(viewZen_norm))
	.subtract(ee.Image(0.034509608392235).multiply(sunZen_norm))
	.add(ee.Image(0.031884395036004).multiply(relAzim_norm));

  return tansig(sum);
}

function neuron5_fapar8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm) {
  var sum =
	ee.Image(1.353023396690570)
	.subtract(ee.Image(0.029929946166941).multiply(b03_norm))
	.add(ee.Image(0.795804414040809).multiply(b04_norm))
	.add(ee.Image(0.348025317624568).multiply(b05_norm))
	.add(ee.Image(0.943567007518504).multiply(b06_norm))
	.subtract(ee.Image(0.276341670431501).multiply(b07_norm))
	.subtract(ee.Image(2.946594180142590).multiply(b8a_norm))
	.add(ee.Image(0.289483073507500).multiply(b11_norm))
	.add(ee.Image(1.044006950440180).multiply(b12_norm))
	.subtract(ee.Image(0.000413031960419).multiply(viewZen_norm))
	.add(ee.Image(0.403331114840215).multiply(sunZen_norm))
	.add(ee.Image(0.068427130526696).multiply(relAzim_norm));

  return tansig(sum);
}

function layer2_fapar8(neuron1_fapar8, neuron2_fapar8, neuron3_fapar8, neuron4_fapar8, neuron5_fapar8) {
  var sum =
	ee.Image(-0.336431283973339)
	.add(ee.Image(2.126038811064490).multiply(neuron1_fapar8))
	.subtract(ee.Image(0.632044932794919).multiply(neuron2_fapar8))
	.add(ee.Image(5.598995787206250).multiply(neuron3_fapar8))
	.add(ee.Image(1.770444140578970).multiply(neuron4_fapar8))
	.subtract(ee.Image(0.267879583604849).multiply(neuron5_fapar8));

  return sum;
}

//---------------------
// FAPAR 3 BAND
//---------------------

exports.get_fapar3band = function(img) {
  img = descale(img, 0.0001)
  var inputoormask = maskinputOOR(img)
  var b03_norm = normalize(img.select('B3'), 0, 0.243425768);
  var b04_norm = normalize(img.select('B4'), 0, 0.297684236);
  var b08_norm = normalize(img.select('B8'), 0.026530282, 0.78139164);
  
  var viewZen_norm = normalize(ee.Image(img.getNumber('MEAN_INCIDENCE_ZENITH_ANGLE_B8')).multiply(degToRad).cos(), 0.918595401, 1);
  var sunZen_norm  = normalize(ee.Image(img.getNumber('MEAN_SOLAR_ZENITH_ANGLE')).multiply(degToRad).cos(), 0.342022871, 0.936206429);
  var relAzim_norm = (ee.Image(img.getNumber('MEAN_SOLAR_AZIMUTH_ANGLE')).subtract(ee.Image(img.getNumber('MEAN_INCIDENCE_AZIMUTH_ANGLE_B8'))).multiply(degToRad)).cos()

  var n1 = neuron1_fapar3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm);
  var n2 = neuron2_fapar3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm);
  var n3 = neuron3_fapar3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm);
  var n4 = neuron4_fapar3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm);
  var n5 = neuron5_fapar3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm);
  
  var l2 = layer2_fapar3(n1, n2, n3, n4, n5);
  
  var fapar = denormalize(l2, 0.000153013, 0.977135097);
  return fapar.updateMask(inputoormask)
}

function neuron1_fapar3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm) {

  var sum =
	ee.Image(-0.019802303)
	.add(ee.Image(1.063928519).multiply(b03_norm))
	.add(ee.Image(0.910752392).multiply(b04_norm))
	.subtract(ee.Image(0.973014301).multiply(b08_norm))
	.subtract(ee.Image(1.26727725).multiply(viewZen_norm))
	.add(ee.Image(0.239696855).multiply(sunZen_norm))
	.subtract(ee.Image(0.837005031).multiply(relAzim_norm));

  return tansig(sum);
  
}

function neuron2_fapar3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm) {
  var sum =
	ee.Image(2.917991233)
	.subtract(ee.Image(1.087124712).multiply(b03_norm))
	.add(ee.Image(2.869208297).multiply(b04_norm))
	.add(ee.Image(0.961199343).multiply(b08_norm))
	.add(ee.Image(0.055681494).multiply(viewZen_norm))
	.subtract(ee.Image(0.267414425).multiply(sunZen_norm))
	.subtract(ee.Image(0.066394844).multiply(relAzim_norm));

  return tansig(sum);
  
}

function neuron3_fapar3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm) {
  var sum =
	ee.Image(- 1.3349831)
	.subtract(ee.Image(0.732287638).multiply(b03_norm))
	.add(ee.Image(0.836483005).multiply(b04_norm))
	.subtract(ee.Image(2.273506421).multiply(b08_norm))
	.add(ee.Image(0.00640356).multiply(viewZen_norm))
	.subtract(ee.Image(0.17567951).multiply(sunZen_norm))
	.subtract(ee.Image(0.022244354).multiply(relAzim_norm));

  return tansig(sum);

}

function neuron4_fapar3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm) {
  
  var sum =
	ee.Image(- 1.38915446)
	.subtract(ee.Image(0.627414923).multiply(b03_norm))
	.add(ee.Image(1.227193715).multiply(b04_norm))
	.subtract(ee.Image(2.532473181).multiply(b08_norm))
	.subtract(ee.Image(0.025617074).multiply(viewZen_norm))
	.subtract(ee.Image(0.125296835).multiply(sunZen_norm))
	.subtract(ee.Image(0.010849463).multiply(relAzim_norm));

  return tansig(sum);

}

function neuron5_fapar3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm) {
  
  var sum =
	ee.Image(0.917074723)
	.add(ee.Image(0.376619209).multiply(b03_norm))
	.add(ee.Image(1.886599724).multiply(b04_norm))
	.subtract(ee.Image(1.841536547).multiply(b08_norm))
	.subtract(ee.Image(0.048726519).multiply(viewZen_norm))
	.add(ee.Image(0.107025026).multiply(sunZen_norm))
	.add(ee.Image(0.005804985).multiply(relAzim_norm));

  return tansig(sum);

}

function layer2_fapar3(neuron1_fapar3, neuron2_fapar3, neuron3_fapar3, neuron4_fapar3, neuron5_fapar3) {

  var sum =
	ee.Image(-0.446230574)
	.add(ee.Image(0.039475758).multiply(neuron1_fapar3))
	.add(ee.Image(0.32828457).multiply(neuron2_fapar3))
	.add(ee.Image(1.149270061).multiply(neuron3_fapar3))
	.subtract(ee.Image(1.610722043).multiply(neuron4_fapar3))
	.subtract(ee.Image(0.733977148).multiply(neuron5_fapar3));

  return sum;
  
}


//---------------------
// LAI 8 BAND
//---------------------

exports.get_lai = function(img) {
  img = descale(img, 0.0001)
  var inputoormask = maskinputOOR(img)
  var b03_norm = normalize(img.select('B3'), 0, 0.253061520471542);
  var b04_norm = normalize(img.select('B4'), 0, 0.290393577911328);
  var b05_norm = normalize(img.select('B5'), 0, 0.305398915248555);
  var b06_norm = normalize(img.select('B6'), 0.006637972542253, 0.608900395797889);
  var b07_norm = normalize(img.select('B7'), 0.013972727018939, 0.753827384322927);
  var b8a_norm = normalize(img.select('B8A'), 0.026690138082061, 0.782011770669178);
  var b11_norm = normalize(img.select('B11'), 0.016388074192258, 0.493761397883092);
  var b12_norm = normalize(img.select('B12'), 0, 0.493025984460231);
  var viewZen_norm = normalize(ee.Image(img.getNumber('MEAN_INCIDENCE_ZENITH_ANGLE_B8')).multiply(degToRad).cos(), 0.918595400582046, 1);
  var sunZen_norm  = normalize(ee.Image(img.getNumber('MEAN_SOLAR_ZENITH_ANGLE')).multiply(degToRad).cos(), 0.342022871159208, 0.936206429175402);
  var relAzim_norm = (ee.Image(img.getNumber('MEAN_SOLAR_AZIMUTH_ANGLE')).subtract(ee.Image(img.getNumber('MEAN_INCIDENCE_AZIMUTH_ANGLE_B8'))).multiply(degToRad)).cos()

  var n1 = neuron1_lai8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm);
  var n2 = neuron2_lai8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm);
  var n3 = neuron3_lai8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm);
  var n4 = neuron4_lai8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm);
  var n5 = neuron5_lai8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm);
  
  var l2 = layer2_lai8(n1, n2, n3, n4, n5);
  
  var lai = denormalize(l2, 0.000319182538301, 14.4675094548151);
  return lai.updateMask(inputoormask)
}


function neuron1_lai8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm) {
  var sum =
	ee.Image(4.96238030555279)
	.subtract(ee.Image(0.023406878966470 ).multiply(b03_norm))
	.add(ee.Image(0.921655164636366 ).multiply(b04_norm))
	.add(ee.Image(0.135576544080099 ).multiply(b05_norm))
	.subtract(ee.Image(1.938331472397950 ).multiply(b06_norm))
	.subtract(ee.Image(3.342495816122680 ).multiply(b07_norm))
	.add(ee.Image(0.902277648009576  ).multiply(b8a_norm))
	.add(ee.Image(0.205363538258614 ).multiply(b11_norm))
	.subtract(ee.Image(0.040607844721716 ).multiply(b12_norm))
	.subtract(ee.Image(0.083196409727092 ).multiply(viewZen_norm))
	.add(ee.Image(0.260029270773809 ).multiply(sunZen_norm))
	.add(ee.Image(0.284761567218845 ).multiply(relAzim_norm));

  return tansig(sum);
  
}

function neuron2_lai8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm) {
  var sum =
	ee.Image(1.416008443981500)
	.subtract(ee.Image(0.132555480856684 ).multiply(b03_norm))
	.subtract(ee.Image(0.139574837333540 ).multiply(b04_norm))
	.subtract(ee.Image(1.014606016898920 ).multiply(b05_norm))
	.subtract(ee.Image(1.330890038649270 ).multiply(b06_norm))
	.add(ee.Image(0.031730624503341 ).multiply(b07_norm))
	.subtract(ee.Image(1.433583541317050 ).multiply(b8a_norm))
	.subtract(ee.Image(0.959637898574699 ).multiply(b11_norm))
	.add(ee.Image(1.133115706551000 ).multiply(b12_norm))
	.add(ee.Image(0.216603876541632 ).multiply(viewZen_norm))
	.add(ee.Image(0.410652303762839 ).multiply(sunZen_norm))
	.add(ee.Image(0.064760155543506 ).multiply(relAzim_norm));
	
  return tansig(sum);
}

function neuron3_lai8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm) {
  var sum =
	ee.Image(1.075897047213310)
	.add(ee.Image(0.086015977724868  ).multiply(b03_norm))
	.add(ee.Image(0.616648776881434 ).multiply(b04_norm))
	.add(ee.Image(0.678003876446556 ).multiply(b05_norm))
	.add(ee.Image(0.141102398644968 ).multiply(b06_norm))
	.subtract(ee.Image(0.096682206883546 ).multiply(b07_norm))
	.subtract(ee.Image(1.128832638862200 ).multiply(b8a_norm))
	.add(ee.Image(0.302189102741375 ).multiply(b11_norm))
	.add(ee.Image(0.434494937299725 ).multiply(b12_norm))
	.subtract(ee.Image(0.021903699490589 ).multiply(viewZen_norm))
  .subtract(ee.Image(0.228492476802263 ).multiply(sunZen_norm))
	.subtract(ee.Image(0.039460537589826 ).multiply(relAzim_norm));
	
  return tansig(sum);
}

function neuron4_lai8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm) {
  var sum =
	ee.Image(1.533988264655420)
	.subtract(ee.Image(0.109366593670404 ).multiply(b03_norm))
	.subtract(ee.Image(0.071046262972729 ).multiply(b04_norm))
	.add(ee.Image(0.064582411478320 ).multiply(b05_norm))
	.add(ee.Image(2.906325236823160 ).multiply(b06_norm))
	.subtract(ee.Image(0.673873108979163 ).multiply(b07_norm))
	.subtract(ee.Image(3.838051868280840 ).multiply(b8a_norm))
	.add(ee.Image(1.695979344531530 ).multiply(b11_norm))
	.add(ee.Image(0.046950296081713 ).multiply(b12_norm))
	.subtract(ee.Image(0.049709652688365 ).multiply(viewZen_norm))
	.add(ee.Image(0.021829545430994 ).multiply(sunZen_norm))
	.add(ee.Image(0.057483827104091 ).multiply(relAzim_norm));

  return tansig(sum);
}

function neuron5_lai8(b03_norm,b04_norm,b05_norm,b06_norm,b07_norm,b8a_norm,b11_norm,b12_norm, viewZen_norm,sunZen_norm,relAzim_norm) {
  var sum =
	ee.Image(3.024115930757230)
	.subtract(ee.Image(0.089939416159969).multiply(b03_norm))
	.add(ee.Image(0.175395483106147).multiply(b04_norm))
	.subtract(ee.Image(0.081847329172620).multiply(b05_norm))
	.add(ee.Image(2.219895367487790).multiply(b06_norm))
	.add(ee.Image(1.713873975136850).multiply(b07_norm))
	.add(ee.Image(0.713069186099534).multiply(b8a_norm))
	.add(ee.Image(0.138970813499201).multiply(b11_norm))
	.subtract(ee.Image(0.060771761518025).multiply(b12_norm))
	.add(ee.Image(0.124263341255473 ).multiply(viewZen_norm))
	.add(ee.Image(0.210086140404351 ).multiply(sunZen_norm))
	.subtract(ee.Image(0.183878138700341 ).multiply(relAzim_norm));

  return tansig(sum);
}

function layer2_lai8(neuron1_lai8, neuron2_lai8, neuron3_lai8, neuron4_lai8, neuron5_lai8) {
  var sum =
	ee.Image(1.096963107077220)
	.subtract(ee.Image(1.500135489728730 ).multiply(neuron1_lai8))
	.subtract(ee.Image(0.096283269121503 ).multiply(neuron2_lai8))
	.subtract(ee.Image(0.194935930577094 ).multiply(neuron3_lai8))
	.subtract(ee.Image(0.352305895755591 ).multiply(neuron4_lai8))
	.add(ee.Image(0.075107415847473 ).multiply(neuron5_lai8));

  return sum;
}


//---------------------
// LAI 3 BAND
//---------------------

exports.get_lai3band = function(img) {
  img = descale(img, 0.0001)
  var inputoormask = maskinputOOR(img)
  var b03_norm = normalize(img.select('B3'), 0, 0.243425768);
  var b04_norm = normalize(img.select('B4'), 0, 0.297684236);
  var b08_norm = normalize(img.select('B8'), 0.026530282, 0.78139164);
  
  var viewZen_norm = normalize(ee.Image(img.getNumber('MEAN_INCIDENCE_ZENITH_ANGLE_B8')).multiply(degToRad).cos(), 0.918595401, 1);
  var sunZen_norm  = normalize(ee.Image(img.getNumber('MEAN_SOLAR_ZENITH_ANGLE')).multiply(degToRad).cos(), 0.342022871, 0.936206429);
  var relAzim_norm = (ee.Image(img.getNumber('MEAN_SOLAR_AZIMUTH_ANGLE')).subtract(ee.Image(img.getNumber('MEAN_INCIDENCE_AZIMUTH_ANGLE_B8'))).multiply(degToRad)).cos()

  var n1 = neuron1_lai3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm);
  var n2 = neuron2_lai3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm);
  var n3 = neuron3_lai3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm);
  var n4 = neuron4_lai3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm);
  var n5 = neuron5_lai3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm);
  
  var l2 = layer2_lai3(n1, n2, n3, n4, n5);
  
  var lai = denormalize(l2, 0.000319183, 14.46750945);
  return lai.updateMask(inputoormask)
}

function neuron1_lai3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm) {

  var sum =
	ee.Image(-1.334407347)
	.add(ee.Image(0.765662655).multiply(b03_norm))
	.add(ee.Image(0.747236156).multiply(b04_norm))
	.subtract(ee.Image(1.441125721).multiply(b08_norm))
	.subtract(ee.Image(0.37266707).multiply(viewZen_norm))
	.subtract(ee.Image(0.805029771).multiply(sunZen_norm))
	.add(ee.Image(1.514114263).multiply(relAzim_norm));

  return tansig(sum);
  
}

function neuron2_lai3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm) {
  var sum =
	ee.Image(-1.259173761)
	.subtract(ee.Image(0.587652462).multiply(b03_norm))
	.subtract(ee.Image(0.809096664).multiply(b04_norm))
	.subtract(ee.Image(0.728566051).multiply(b08_norm))
	.subtract(ee.Image(0.803710774).multiply(viewZen_norm))
	.add(ee.Image(0.093313347).multiply(sunZen_norm))
	.add(ee.Image(1.347624019).multiply(relAzim_norm));

  return tansig(sum);
  
}

function neuron3_lai3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm) {
  var sum =
	ee.Image(0.734655902)
	.subtract(ee.Image(0.156091129).multiply(b03_norm))
	.add(ee.Image(1.426190194).multiply(b04_norm))
	.subtract(ee.Image(0.840785867).multiply(b08_norm))
	.add(ee.Image(0.067236054).multiply(viewZen_norm))
	.subtract(ee.Image(0.32100852).multiply(sunZen_norm))
	.subtract(ee.Image(0.190242198).multiply(relAzim_norm));

  return tansig(sum);

}

function neuron4_lai3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm) {
  var sum =
	ee.Image(1.869059826)
	.subtract(ee.Image(0.214752449).multiply(b03_norm))
	.subtract(ee.Image(0.067188443).multiply(b04_norm))
	.subtract(ee.Image(1.636409461).multiply(b08_norm))
	.add(ee.Image(0.028918314).multiply(viewZen_norm))
	.add(ee.Image(0.126143729).multiply(sunZen_norm))
	.add(ee.Image(0.104116201).multiply(relAzim_norm));

  return tansig(sum);

}

function neuron5_lai3(b03_norm,b04_norm,b08_norm,viewZen_norm,sunZen_norm,relAzim_norm) {
  var sum =
	ee.Image(1.594139175)
	.add(ee.Image(1.715027862).multiply(b03_norm))
	.add(ee.Image(0.871304754).multiply(b04_norm))
	.subtract(ee.Image(2.027147985).multiply(b08_norm))
	.subtract(ee.Image(0.055184617).multiply(viewZen_norm))
	.add(ee.Image(0.008336523).multiply(sunZen_norm))
	.add(ee.Image(0.006279199).multiply(relAzim_norm));

  return tansig(sum);

}

function layer2_lai3(neuron1_lai3, neuron2_lai3, neuron3_lai3, neuron4_lai3, neuron5_lai3) {
  var sum =
	ee.Image(0.480722333)
	.add(ee.Image(0.017013185).multiply(neuron1_lai3))
	.subtract(ee.Image(0.044391605).multiply(neuron2_lai3))
	.subtract(ee.Image(0.186661848).multiply(neuron3_lai3))
	.subtract(ee.Image(1.090955074).multiply(neuron4_lai3))
	.subtract(ee.Image(0.260050915).multiply(neuron5_lai3));

  return sum;
  
}



//------------------
// GENERAL FUNCTIONS
//------------------

function descale(scaled, scalefactor) {
  return ee.Image(scaled.multiply(scalefactor).copyProperties(scaled))
}

function normalize(unnormalized, min, max) {
  return ee.Image(2).multiply(unnormalized.subtract(ee.Image(min))).divide(
          ee.Image(max).subtract(ee.Image(min))).subtract(ee.Image(1))
          .copyProperties(unnormalized)
}
function denormalize(normalized, min, max) {
  return ee.Image(0.5).multiply(normalized.add(ee.Image(1))).multiply(ee.Image(max).subtract(ee.Image(min))).add(ee.Image(min))
}

function tansig(input) {
  return ee.Image(2).divide(ee.Image(1).add((ee.Image(-2).multiply(input)).exp())).subtract(ee.Image(1))
  .copyProperties(input)
}

function maskinputOOR(img) {
  var mask = img.select('B3').lte(0.26).and(
    img.select('B4').lte(0.30)).and(
      img.select('B5').lte(0.32)).and(
        img.select('B6').lte(0.62)).and(
          img.select('B7').lte(0.75)).and(
            img.select('B8A').lte(0.78)).and(
              img.select('B11').lte(0.52)).and(
                img.select('B12').lte(0.50));
    
  return mask
    
}
