/* 
 * The MIT License
 *
 * Copyright 2018 Bruce Schubert.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

/* global expect */

import Rothermel from '../src/Rothermel';

const TO_RADIANS = Math.PI / 180;
const TO_DEGREES = 180 / Math.PI;

function toRadians(degrees) {
  return degrees * TO_RADIANS;
}

describe('Rothermel', () => {
  it('should exist', () => {
    expect(Rothermel).toBeTruthy();
  });
});


describe('meanBulkDensity', () => {
  const w0 = [1, 2, 3];  // sum = 6
  const height = 2;
  const rho_b = 3;

  // expect().toThrowError requires the function to be tested be wrapped in an anonymous function
  it('throws with height zero', () => {
    expect(() => {
      Rothermel.meanBulkDensity(w0, 0);
    }).toThrowError(RangeError);
  });

  it('throws with height < zero', () => {
    expect(() => {
      Rothermel.meanBulkDensity(w0, -1);
    }).toThrowError(RangeError);
  });

  it('result is valid', () => {
    expect(Rothermel.meanBulkDensity(w0, height)).toBe(rho_b);
  });
});


describe('meanPackingRatio', () => {
  const rho_b = 1;
  const rho_p = 10;
  const beta = 0.1;

  // Note: expect().toThrowError requires the function to be tested be wrapped in an anonymous function

  it('throws with rho_p zero', () => {
    expect(() => {
      Rothermel.meanPackingRatio(rho_b, 0);
    }).toThrowError(RangeError);
  });

  it('throws with rho_p < zero', () => {
    expect(() => {
      Rothermel.meanPackingRatio(rho_b, -1);
    }).toThrowError(RangeError);
  });

  it('throws with beta < zero', () => {
    expect(() => {
      Rothermel.meanPackingRatio(rho_b, -1);
    }).toThrowError(RangeError);
  });

  it('throws with beta > 0.12', () => {
    expect(() => {
      Rothermel.meanPackingRatio(rho_b, 7.6);
    }).toThrowError(RangeError);
  });

  it('throws with beta < 0', () => {
    expect(() => {
      Rothermel.meanPackingRatio(-1, rho_p);
    }).toThrowError(RangeError);
  });

  it('result is valid', () => {
    expect(Rothermel.meanPackingRatio(rho_b, rho_p)).toBe(beta);
  });
});


describe('optimalPackingRatio', () => {
  const sigma = 1;
  const beta_opt = 3.348 * Math.pow(sigma, -0.8189);

  it('algorithm has not changed', () => {
    expect(Rothermel.optimalPackingRatio(sigma)).toBeCloseTo(beta_opt, 5);
  });
});


describe('characteristicSAV', () => {
  const sv = [3];  // ft2/ft3
  const w0 = [2];  // lbs/ft2
  const sigma = (sv[0] * sv[0] * w0[0]) / (sv[0] * w0[0]);

  it('algorithm has not changed', () => {
    expect(Rothermel.characteristicSAV(sv, w0)).toBe(sigma);
  });
});


describe('reactionVelocity', () => {
  const sigma = 0.3;      // SAV ft2/ft3
  const beta_ratio = 1; // packing ratio
  const gamma = 0.00033194651982916144;

  it('algorithm has not changed', () => {
    expect(Rothermel.reactionVelocity(sigma, beta_ratio)).toBeCloseTo(gamma, 10);
  });
});


describe('reactionIntensity', () => {
  const gamma = 2;
  const heat = 3;
  const eta_M = 4;
  const eta_s = 5;
  const I_r = gamma * heat * eta_M * eta_s;

  it('algorithm has not changed', () => {
    expect(Rothermel.reactionIntensity(gamma, heat, eta_M, eta_s)).toBe(I_r);
  });
});


describe('flameResidenceTime', () => {
  const sigma = 2;
  const tau = 384 / sigma;

  it('algorithm has not changed', () => {
    expect(Rothermel.flameResidenceTime(sigma)).toBe(tau);
  });
});


describe('heatRelease', () => {
  const I_r = 2;
  const tau = 3;
  const hpa = I_r * tau;

  it('algorithm has not changed', () => {
    expect(Rothermel.heatRelease(I_r, tau)).toBe(hpa);
  });
});


describe('propagatingFluxRatio', () => {
  const sigma = 2;
  const beta = 3;
  const xi = Math.exp((0.792 + 0.681 * Math.sqrt(sigma)) * (beta + 0.1)) / (192 + 0.2595 * sigma);

  it('algorithm has not changed', () => {
    expect(Rothermel.propagatingFluxRatio(sigma, beta)).toBeCloseTo(xi, 10);
  });
});


describe('effectiveHeatingNumber', () => {
  const sv = 2;
  const epsilon = Math.exp(-138. / sv);

  it('algorithm has not changed', () => {
    expect(Rothermel.effectiveHeatingNumber(sv)).toBe(epsilon);
  });
});


describe('heatOfPreignition', () => {
  const Mf = 2;
  const Q_ig = 250 + 1116 * (Mf * 0.01); // Mf = [fraction]

  it('algorithm has not changed', () => {
    expect(Rothermel.heatOfPreignition(Mf)).toBe(Q_ig);
  });
});


describe('heatSink', () => {
  const preignitionHeat = [2];
  const effectiveHeating = [3];
  const sw = [4];
  const density = 5;
  const hsk = density * ((preignitionHeat[0] * effectiveHeating[0] * sw[0]) / sw[0]);

  it('algorithm has not changed', () => {
    expect(Rothermel.heatSink(preignitionHeat, effectiveHeating, sw, density)).toBe(hsk);
  });
});


describe('windFactor', () => {
  const midFlameWindSpd = 2;
  const sigma = 0.3;
  const beta_ratio = 0.4;
  const phi_w = 13.551811735940076;

  it('algorithm has not changed', () => {
    expect(Rothermel.windFactor(midFlameWindSpd, sigma, beta_ratio)).toBeCloseTo(phi_w, 10);
  });
});


describe('slopeFactor', () => {
  const slopeDegrees = 45;
  const beta = 0.3;
  const phi_s = 7.569829322728009;

  it('algorithm has not changed', () => {
    expect(Rothermel.slopeFactor(slopeDegrees, beta)).toBeCloseTo(phi_s, 10);
  });
});


describe('effectiveWindSpeed', () => {
  const phiEw = 10;
  const beta_ratio = 0.6;
  const sigma = 0.7;
  const efw = 5.58851070023703; //[ft/min]

  it('algorithm has not changed', () => {
    expect(Rothermel.effectiveWindSpeed(phiEw, beta_ratio, sigma)).toBeCloseTo(efw, 10);
  });
});


describe('rateOfSpread', () => {
  const reactionIntensity = 2;
  const propogatingFlux = 3;
  const windFactor = 4;
  const slopeFactor = 5;
  const heatSink = 6;
  const ros = (reactionIntensity * propogatingFlux * (1 + windFactor + slopeFactor)) / heatSink; // 10 [ft/min]

  it('algorithm has not changed', () => {
    expect(Rothermel.rateOfSpread(reactionIntensity, propogatingFlux, windFactor, slopeFactor, heatSink)).toBe(ros);
  });
});


describe('rateOfSpreadNoWindNoSlope', () => {
  const reactionIntensity = 2;
  const propogatingFlux = 3;
  const heatSink = 6;
  const ros = (reactionIntensity * propogatingFlux) / heatSink; // 1 [ft/min]

  it('algorithm has not changed', () => {
    expect(Rothermel.rateOfSpreadNoWindNoSlope(reactionIntensity, propogatingFlux, heatSink)).toBe(ros);
  });
});


describe('flameZoneDepth', () => {
  const rateOfSpread = 2;
  const flameResidenceTime = 3;
  const fzd = rateOfSpread * flameResidenceTime;

  it('algorithm has not changed', () => {
    expect(Rothermel.flameZoneDepth(rateOfSpread, flameResidenceTime)).toBe(fzd);
  });
});


describe('firelineIntensity', () => {
  const reactionIntensity = 3;
  const flameZoneDepth = 7;
  const I = reactionIntensity * flameZoneDepth / 60;

  it('algorithm has not changed', () => {
    expect(Rothermel.firelineIntensity(flameZoneDepth, reactionIntensity)).toBeCloseTo(I, 10);
  });
});


describe('flameLength', () => {
  const firelineIntensity = 100;
  const L = 3.74293696996202;

  it('algorithm has not changed', () => {
    expect(Rothermel.flameLength(firelineIntensity)).toBeCloseTo(L, 10);
  });
});


describe('midFlameWindAdjustmentFactor', () => {
  const fuelDepth = 7;
  const waf = 0.5703218562580741;

  it('algorithm has not changed', () => {
    expect(Rothermel.midFlameWindAdjustmentFactor(fuelDepth)).toBeCloseTo(waf, 10);
  });
});


describe('calcWindSpeedMidFlame', () => {
  const wndSpd20Ft = 10;
  const fuelDepth = 7;
  const waf = 5.703218562580741;

  it('algorithm has not changed', () => {
    expect(Rothermel.calcWindSpeedMidFlame(wndSpd20Ft, fuelDepth)).toBeCloseTo(waf, 10);
  });
});


describe('calcWindSpeedNearFuel', () => {
  const wndSpd20Ft = 10;
  const fuelDepth = 7;
  const U_h = 3.1165128757271803;

  it('algorithm has not changed', () => {
    expect(Rothermel.calcWindSpeedNearFuel(wndSpd20Ft, fuelDepth)).toBeCloseTo(U_h, 10);
  });
});


describe('eccentricity', () => {
  const effectiveWind = 10;
  const eccentricity = 0.23342132431186122;

  it('algorithm has not changed', () => {
    expect(Rothermel.eccentricity(effectiveWind)).toBeCloseTo(eccentricity, 10);
  });
});


describe('calcFuelTemp', () => {
  const I = 3;
  const T_a = 67; // F
  const U_h = 7; // mph
  const T_f = 89.90076335877862;

  it('algorithm has not changed', () => {
    expect(Rothermel.calcFuelTemp(I, T_a, U_h)).toBeCloseTo(T_f, 10);
  });
});


describe('calcRelativeHumidityNearFuel', () => {
  const H_a = 50; // %
  const T_a = 67; // F
  const T_f = 89; // F
  const H_f = 24.19202432339955; // %

  it('algorithm has not changed', () => {
    expect(Rothermel.calcRelativeHumidityNearFuel(H_a, T_f, T_a)).toBeCloseTo(H_f, 10);
  });
});


describe('calcEarthSunDistanceSqrd', () => {

  it('algorithm has not changed', () => {
    const delta = 34 * TO_RADIANS;
    const r2 = 1.047651;
    expect(Rothermel.calcEarthSunDistanceSqrd(delta)).toBeCloseTo(r2, 5);
  });

  it('matches July 5 published result', () => {
    const NJ = Rothermel.calcJulianDate(7, 5, 2000);
    const delta = Rothermel.calcSolarDeclinationAngle(NJ);
    const r = 1.01671;        // from text for July 5´
    const expResult = r * r;  // r-squared
    
    expect(Rothermel.calcEarthSunDistanceSqrd(delta)).toBeCloseTo(expResult, 2);
  });
  it('matches Jan 3 expected result', () => {
    const NJ = Rothermel.calcJulianDate(1, 3, 2000);
    const delta = Rothermel.calcSolarDeclinationAngle(NJ);
    const r = 0.98324;        // from text for Jan 3
    const expResult = r * r;  // r-squared
    
    expect(Rothermel.calcEarthSunDistanceSqrd(delta)).toBeCloseTo(expResult, 2);
  });
});


describe('calcSolarIrradianceOnHorzSurface', () => {
  const I_a = 3;
  const r2 = 1.047651;
  const A = 45 * TO_RADIANS;
  const I = 2.0248349341141676;

  it('algorithm has not changed', () => {
    expect(Rothermel.calcSolarIrradianceOnHorzSurface(I_a, r2, A)).toBeCloseTo(I, 10);
  });
});


describe('calcOpticalAirMass', () => {

  it('zenith at sea level', () => {
    const A = toRadians(90);  // zenith
    const E = 0.0;            // sea level [feet]
    const M = 1;              // secant A
    expect(Rothermel.calcOpticalAirMass(A, E)).toBeCloseTo(M, 3);
  });
  
  it('decrease solar angle', () => {
    // test solar angle, as angle goes down, ratio goes up
    const A = toRadians(45);
    const E = 0.0;            // sea level [feet]
    const M = 1.414;          // 1/sin(A)
    expect(Rothermel.calcOpticalAirMass(A, E)).toBeCloseTo(M, 3);
  });
  
  it('increase elevation', () => {
    // test elevation, as elevation increases, ratio goes down
    const A = toRadians(45); 
    const E = 3280 * 5;     // 5km [feet]
    const M = 0.707;        // @ 5km, you are above ~1/2 the air mass
    expect(Rothermel.calcOpticalAirMass(A, E)).toBeCloseTo(M, 1);
  });
});


describe('calcAttenuatedIrradiance', () => {
  const M = 1.3;
  const S_c = 50;
  const p = 1; // ft
  const I_a = 0.99;

  it('algorithm has not changed', () => {
    expect(Rothermel.calcAttenuatedIrradiance(M, S_c, p)).toBeCloseTo(I_a, 3);
  });
});


describe('calcIrradianceOnASlope', () => {
  const alpha = 45 * TO_RADIANS;  // slope angle
  const beta = 135 * TO_RADIANS;  // aspect
  const A = 60 * TO_RADIANS;      // solar altitude
  const Z = 180 * TO_RADIANS;     // solar azimuth
  const I_a = 0.99;               // attenuated irradiance
  const I = 0.8537487113388367;

  it('algorithm has not changed', () => {
    expect(Rothermel.calcIrradianceOnASlope(alpha, beta, A, Z, I_a)).toBeCloseTo(I, 10);
  });
});


describe('calcCanadianStandardDailyFineFuelMoisture', () => {
  const m_0 = 20;   // initial moisture at noon
  const T_f = 67;   // air temp next to fuel
  const H_f = 50;   // RH next to fuel
  const W = 10;     // 20ft winds mph
  const R = 0;      // rainfall
  const fm = 15.550628477777991;

  it('algorithm has not changed', () => {
    expect(Rothermel.calcCanadianStandardDailyFineFuelMoisture(m_0, T_f, H_f, W, R)).toBeCloseTo(fm, 10);
  });
});


describe('calcCanadianHourlyFineFuelMoisture', () => {
  const m_0 = 20;   // previous hour fuel moisture
  const T_c = 20;   // air temp C
  const H = 50;     // RH 
  const W_k = 16;   // 20ft wind speed KPH
  const fm = 19.07905001378493;

  it('algorithm has not changed', () => {
    expect(Rothermel.calcCanadianHourlyFineFuelMoisture(m_0, H, T_c, W_k)).toBeCloseTo(fm, 10);
  });
});


describe('calcFineDeadFuelMoisture', () => {
  const m_0 = 20;   // initial fuel moisture
  const T_c = 20;   // air temp C
  const H = 50;     // RH 
  const fm = 13.68836723429169;

  it('algorithm has not changed', () => {
    expect(Rothermel.calcFineDeadFuelMoisture(m_0, T_c, H)).toBeCloseTo(fm, 10);
  });
});


describe('calcJulianDate', () => {
  it('1/32/2000 should be 32', () => {
    expect(Rothermel.calcJulianDate(1, 32, 2000)).toBe(32);
  });
  it('12/31/2000 should be 366', () => {
    expect(Rothermel.calcJulianDate(12, 31, 2000)).toBe(366);
  });
  it('12/31/2009 should be 365', () => {
    expect(Rothermel.calcJulianDate(12, 31, 2009)).toBe(365);
  });
});


describe('calcLocalHourAngle', () => {
  it('morning', () => {
    const localTime = 6.0;
    const h = 0;
    expect(Rothermel.calcLocalHourAngle(localTime)).toBeCloseTo(h, 5);
  });
  it('local noon', () => {
    const localTime = 12.0;
    const h = (90 * TO_RADIANS);
    expect(Rothermel.calcLocalHourAngle(localTime)).toBeCloseTo(h, 5);
  });
  it('evening', () => {
    const localTime = 18.0;
    const h = (180 * TO_RADIANS);
    expect(Rothermel.calcLocalHourAngle(localTime)).toBeCloseTo(h, 5);
  });
  it('midnight', () => {
    const localTime = 24.0;
    const h = (270 * TO_RADIANS);
    expect(Rothermel.calcLocalHourAngle(localTime)).toBeCloseTo(h, 5);
  });
});


describe('calcSolarDeclinationAngle', () => {

  it('winter solstice should be -23.5', () => {
    const NJ = Rothermel.calcJulianDate(12, 21, 2009);
    const expResult = (-23.5 * TO_RADIANS);
    expect(Rothermel.calcSolarDeclinationAngle(NJ)).toBeCloseTo(expResult, 5);
  });

  it('summer solstice should be 23.5', () => {
    const NJ = Rothermel.calcJulianDate(6, 21, 2009);
    const expResult = (23.5 * TO_RADIANS);
    expect(Rothermel.calcSolarDeclinationAngle(NJ)).toBeCloseTo(expResult, 5);
  });

  it('spring equinox should be 0', () => {
    const NJ = Rothermel.calcJulianDate(3, 22, 2009);
    const expResult = 0;
    expect(Rothermel.calcSolarDeclinationAngle(NJ)).toBeCloseTo(expResult, 2);
  });

  it('fall equinox should be 0', () => {
    const NJ = Rothermel.calcJulianDate(9, 21, 2009);
    const expResult = 0;
    expect(Rothermel.calcSolarDeclinationAngle(NJ)).toBeCloseTo(expResult, 2);
  });

});


describe('calcSolarAltitudeAngle', () => {
  const phi = toRadians(0);  // equator
  const delta = Rothermel.calcSolarDeclinationAngle(Rothermel.calcJulianDate(3, 21, 2000));
  it('3/21/2000 @ 0600', () => {
    const h = Rothermel.calcLocalHourAngle(6.0);   // local time
    const A = toRadians(0);
    expect(Rothermel.calcSolarAltitudeAngle(h, phi, delta)).toBeCloseTo(A, 5);
  });
  it('3/21/2000 @ 1200', () => {
    const h = Rothermel.calcLocalHourAngle(12.0);   // local time
    const A = toRadians(90);
    expect(Rothermel.calcSolarAltitudeAngle(h, phi, delta)).toBeCloseTo(A, 5);
  });
  it('3/21/2000 @ 1800', () => {
    const h = Rothermel.calcLocalHourAngle(18.0);   // local time
    const A = toRadians(0);
    expect(Rothermel.calcSolarAltitudeAngle(h, phi, delta)).toBeCloseTo(A, 5);
  });
});



describe('calcSolarAzimuthAngle', () => {
  const phi = toRadians(-34.2);                         // ventura
  const NJ = Rothermel.calcJulianDate(3, 21, 2000);     // vernal equinox
  const delta = Rothermel.calcSolarDeclinationAngle(NJ);
  
  it('3/21/2000 @ 0600', () => {
    const h = Rothermel.calcLocalHourAngle(6.0);        // morning - local time
    const A = Rothermel.calcSolarAltitudeAngle(h, phi, delta);
    const Z  = toRadians(360);
    expect(Rothermel.calcSolarAzimuthAngle(h, phi, delta, A)).toBeCloseTo(Z, 5);
  });
  it('3/21/2000 @ 1200', () => {
    const h = Rothermel.calcLocalHourAngle(12.0);        // noon - local time
    const A = Rothermel.calcSolarAltitudeAngle(h, phi, delta);
    const Z  = toRadians(90);
    expect(Rothermel.calcSolarAzimuthAngle(h, phi, delta, A)).toBeCloseTo(Z, 5);
  });
  it('3/21/2000 @ 1800', () => {
    const h = Rothermel.calcLocalHourAngle(18.0);        // noon - local time
    const A = Rothermel.calcSolarAltitudeAngle(h, phi, delta);
    const Z  = toRadians(180);
    expect(Rothermel.calcSolarAzimuthAngle(h, phi, delta, A)).toBeCloseTo(Z, 5);
  });

});

