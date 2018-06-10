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

import Rothermel from '../src/Rothermel'

const TO_RADIANS = Math.PI / 180;
const TO_DEGREES = 180 / Math.PI;

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
    expect(Rothermel.optimalPackingRatio(sigma)).toBe(beta_opt);
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
  const sigma = 2;
  const beta_ratio = 3;
  const gamma = 4.702292202212657e-33;

  it('algorithm has not changed', () => {
    expect(Rothermel.reactionVelocity(sigma, beta_ratio)).toBe(gamma);
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
    expect(Rothermel.propagatingFluxRatio(sigma, beta)).toBe(xi);
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
    expect(Rothermel.windFactor(midFlameWindSpd, sigma, beta_ratio)).toBe(13.551811735940076);
  });
});


describe('slopeFactor', () => {
  const slopeDegrees = 45;
  const beta = 0.3;
  const phi_s = 7.569829322728009;
  
  it('algorithm has not changed', () => {
    expect(Rothermel.slopeFactor(slopeDegrees, beta)).toBe(phi_s);
  });
});


describe('effectiveWindSpeed', () => {
  const phiEw = 10;
  const beta_ratio = 0.6;
  const sigma = 0.7;
  const efw = 5.58851070023703; //[ft/min]
  
  it('algorithm has not changed', () => {
    expect(Rothermel.effectiveWindSpeed(phiEw, beta_ratio, sigma)).toBe(efw);
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
    expect(Rothermel.firelineIntensity(flameZoneDepth, reactionIntensity)).toBe(I);
  });
});


describe('flameLength', () => {
  const firelineIntensity = 100;
  const L = 3.74293696996202;
  
  it('algorithm has not changed', () => {
    expect(Rothermel.flameLength(firelineIntensity)).toBe(L);
  });
});


describe('midFlameWindAdjustmentFactor', () => {
  const fuelDepth = 7;
  const waf = 0.5703218562580741;
  
  it('algorithm has not changed', () => {
    expect(Rothermel.midFlameWindAdjustmentFactor(fuelDepth)).toBe(waf);
  });
});


describe('calcWindSpeedMidFlame', () => {
  const wndSpd20Ft = 10;
  const fuelDepth = 7;
  const waf = 5.703218562580741;
  
  it('algorithm has not changed', () => {
    expect(Rothermel.calcWindSpeedMidFlame(wndSpd20Ft, fuelDepth)).toBe(waf);
  });
});


describe('calcWindSpeedNearFuel', () => {
  const wndSpd20Ft = 10;
  const fuelDepth = 7;
  const U_h = 3.1165128757271803;
  
  it('algorithm has not changed', () => {
    expect(Rothermel.calcWindSpeedNearFuel(wndSpd20Ft, fuelDepth)).toBe(U_h);
  });
});


describe('eccentricity', () => {
  const effectiveWind = 10;
  const eccentricity = 0.23342132431186122;
  
  it('algorithm has not changed', () => {
    expect(Rothermel.eccentricity(effectiveWind)).toBe(eccentricity);
  });
});


describe('calcFuelTemp', () => {
  const I = 3;
  const T_a = 67; // F
  const U_h = 7; // mph
  const T_f = 89.90076335877862;
  
  it('algorithm has not changed', () => {
    expect(Rothermel.calcFuelTemp(I, T_a, U_h) ).toBe(T_f);
  });
});


describe('calcRelativeHumidityNearFuel', () => {
  const H_a = 50; // %
  const T_a = 67; // F
  const T_f = 89; // F
  const H_f = 24.19202432339955; // %
  
  it('algorithm has not changed', () => {
    expect(Rothermel.calcRelativeHumidityNearFuel(H_a, T_f, T_a)).toBe(H_f);
  });
});


describe('calcEarthSunDistanceSqrd', () => {
  const delta = 34 * TO_RADIANS;
  const r2 = 1.047651;
  
  it('algorithm has not changed', () => {
    expect(Rothermel.calcEarthSunDistanceSqrd(delta)).toBe(r2);
  });
});


describe('calcSolarIrradianceOnHorzSurface', () => {
  const I_a = 3;
  const r2 = 1.047651;
  const A = 45 * TO_RADIANS;
  const I = 2.0248349341141676; 
  
  it('algorithm has not changed', () => {
    expect(Rothermel.calcSolarIrradianceOnHorzSurface(I_a, r2, A)).toBe(I);
  });
});


describe('calcOpticalAirMass', () => {
  const A = 45 * TO_RADIANS;
  const E = 1000; // ft
  const M = 1.3522550283652721; 
  
  it('algorithm has not changed', () => {
    expect(Rothermel.calcOpticalAirMass(A, E)).toBe(M);
  });
});


describe('calcAttenuatedIrradiance', () => {
  const M = 1.3;
  const S_c = 50;
  const p = 1; // ft
  const I_a = 0.99; 
  
  it('algorithm has not changed', () => {
    expect(Rothermel.calcAttenuatedIrradiance(M, S_c, p) ).toBe(I_a);
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
    expect(Rothermel.calcIrradianceOnASlope(alpha, beta, A, Z, I_a)).toBe(I);
  });
});