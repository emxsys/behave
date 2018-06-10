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

const HALF_PI = Math.PI / 2;
const TO_RADIANS = Math.PI / 180;
const TO_DEGREES = 180 / Math.PI;

/**
 * The Rothermel, et al, fire spread model developed for modeling wildland fire behavior.
 *
 * References:
 * <ul>
 * <li><a name="bib_1000"></a>Albini, F.A., 1976, Estimating Wildfire Behavior and Effects, General
 * Technical Report INT-30, USDA Forest Service, Intermountain Forest and Range Experiment Station
 *
 * <li><a name="bib_1001"></a>Anderson, H.A., 1983, Predicting Wind-driven Wild Land Fire Size and
 * Shape, Research Paper INT-305, USDA Forest Service, Intermountain Forest and Range Experiment
 * Station
 *
 * <li><a name="bib_1002"></a>Anderson, K., 2009, A Comparison of Hourly Fire Fuel Moisture Code
 * Calculations within Canada, Canadian Forest Service
 *
 * <li><a name="bib_1003"></a>Rothermel, R.C., 1972, A mathematical model for predicting fire spread
 * in wildland fuels, General Technical Report INT-115, USDA Forest Service, Intermountain Forest
 * and Range Experiment Station
 *
 * <li><a name="bib_1004"></a>Rothermel et al, 1986, Modeling Moisture Content of Fine Dead Wildland
 * Fuels: Input to the BEHAVE Fire Prediction System, Research Paper, INT-359, USDA Forest Service,
 * Intermountain Research Station
 *
 * <li><a name="bib_1005"></a>Van Wagner, C.E., 1977, A Method of Computing Fine Fuel Moisture
 * Content Throughout the Diurnal Cycle, Information Report PS-X-69, Petawawa Forest Experiment
 * Station, Canadian Forest Service
 * </ul>
 * Other Sources:
 * <ul>
 * <li>BehavePlus5, xfblib.cpp, Copyright Collin D. Bevins.
 * <li>Firelib v1.04, firelib.c, Copyright Collin D. Bevins.
 * </ul>
 *
 * @author Bruce Schubert
 */
export default class Rothermel {

  /**
   * Calculates the mean bulk density (fuel-bed weight per unit volume): rho_b.
   * 
   * @param {Number[]} w0 An array of fuel particle loading values [lb/ft2]
   * @param {Number} height The fuel bed height [ft]
   * 
   * @return {Number} rho_b [lbs/ft3]
   */
  static meanBulkDensity(w0, height) {
    if (height <= 0) {
      throw new RangeError("height must be > 0.")
    }
    const w0_t = w0.reduce((accumulator, currentValue) => accumulator + currentValue)
    const rho_b = w0_t / height

    return rho_b
  }

  /**
   * Calculates the mean packing ratio for the fuel: beta.
   *
   * The compactness of the fuel bed is defined by the packing ratio, which is defined as the
   * fraction of the fuel array volume that is occupied by the fuel. Rothermel 1972: eq. (74)
   *
   * @param {Number} rho_b The mean bulk density of the fuel bed [lbs/ft3].
   * @param {Number} rho_p The oven-dry fuel-particle density [lbs/ft3].
   *
   * @return {Number} beta [dimensionless]
   */
  static meanPackingRatio(rho_b, rho_p) {
    if (rho_p <= 0) {
      throw new RangeError("rho_p must be > 0.");
    }
    let beta = rho_b / rho_p;
    if (beta > 0.12 || beta < 0) {
      throw new RangeError("Mean packing ratio [beta] out of limits [0,0.12]: " + beta);
    }
    return beta;
  }

  /**
   * Computes the optimal packing ratio for the fuel: beta_opt.
   *
   * Optimum packing ratio is a term used in the Rothermel's (1972) surface fire spread model
   * indicating the packing ratio that optimizes the reaction velocity term of the spread model.
   * Optimum packing ratio is a function of the fineness of fuel particles, which is measured by
   * the characteristic surface-area-to-volume ratio of the fuelbed. Optimum packing ratio does
   * not optimize fire behavior (rate of spread or fireline intensity). Fire Science Glossary
   * [electronic]. http://www.firewords.net
   *
   * @param {Number} sigma The characteristic SAV ratio for the fuel complex [ft2/ft3].
   *
   * @return {Number} beta_opt [dimensionless]
   */
  static optimalPackingRatio(sigma) {
    let beta_opt = 3.348 * Math.pow(sigma, -0.8189)
    return beta_opt
  }

  /**
   * Computes the characteristic surface-area-to-volume ratio for the fuel complex: sigma.
   *
   * In Rothermel's (1972) surface fire spread model, characteristic surface-area-to-volume (SAV)
   * ratio constitutes the fuelbed-average SAV weighted by particle surface area. Surface-area
   * weighting emphasizes fine fuel because finer fuel particles have larger SAV ratios. Fire
   * Science Glossary [electronic]. http://www.firewords.net
   *
   * Rothermel 1972: eq. (71) and (72).
   *
   * @param {Number[]} sv An array of fuel particle SAV ratio values [ft2/ft3].
   * @param {Number[]} w0 An array of fuel particle loading values [lbs/ft2].
   *
   * @return {Number} sigma [ft2/ft3]
   */
  static characteristicSAV(sv, w0) {
    let sw_t = 0      // sw = (sv * w)
    let s2w_t = 0     // s2w = (sv^2 * w)
    const len = sv.length
    for (let i = 0; i < len; i++) {
      sw_t += sv[i] * w0[i]
      s2w_t += sv[i] * sv[i] * w0[i]
    }
    if (sw_t <= 0) {
      throw new RangeError("w0 total loading must be > 0.");
    }
    const sigma = s2w_t / sw_t
    return sigma;
  }

  /**
   * Calculates the potential reaction velocity: gamma.
   *
   * Rothermel 1972: eq. (68),(70) and Albini 1976: pg. 88
   *
   * @param {Number} sigma The characteristic SAV ratio [ft2/ft3].
   * @param {Number} beta_ratio The relative packing ratio [beta/beta_opt].
   * @return {Number} gamma [1/min]
   */
  static reactionVelocity(sigma, beta_ratio) {
    const sigma15 = Math.pow(sigma, 1.5)
    const A = 133 / Math.pow(sigma, 0.7913)    // Albini 
    const gamma_max = sigma15 / (495 + 0.0594 * sigma15)
    const gamma = gamma_max * Math.pow(beta_ratio, A) * Math.exp(A * (1 - beta_ratio))
    return gamma
  }

  /**
   * Calculates the reaction intensity: I_r.
   *
   * The rate of heat release, per unit area of the flaming fire front, expressed as heat
   * energy/area/time, such as Btu/square foot/minute, or Kcal/square meter/second.
   *
   * Rothermel 1972: eq. (58), (59) thru (60)
   *
   * @param {Number} gamma The potential reaction velocity.
   * @param {Number} heat The low heat content [Btu/lb]
   * @param {Number} eta_M The moisture damping coefficient.
   * @param {Number} eta_s The mineral damping coefficient.
   *
   * @return {Number} I_r [BTU/ft2/min].
   */
  static reactionIntensity(gamma, heat, eta_M, eta_s) {
    const I_r = gamma * heat * eta_M * eta_s;
    return I_r;
  }

  /**
   * Calculates the flame residence time: tau.
   *
   * Albini (1976): p.91
   *
   * @param {Number} sigma The characteristic SAV ratio [ft2/ft3].
   * @return {Number} tau [min].
   */
  static flameResidenceTime(sigma) {
    if (sigma <= 0) {
      throw new RangeError("sigma must be > 0.")
    }
    const tau = 384 / sigma
    return tau
  }

  /**
   * Calculates the heat release per unit area: hpa.
   *
   * @param {Number} I_r The reaction intensity [Btu/ft2/min].
   * @param {Number} tau The flame residence time [min].
   *
   * @return {NUmber} hpa [Btu/ft2]
   */
  static heatRelease(I_r, tau) {
    const hpa = I_r * tau
    return hpa
  }

  /**
   * Gets the propagating flux ratio: xi.
   *
   * The no-wind propagating flux ratio is a function of the mean packing ratio (beta) and the
   * characteristic SAV ratio (sigma).
   *
   * Rothermel 1972: eq. (42)(76)
   *
   * @param {Number} sigma The characteristic SAV ratio [ft2/ft3].
   * @param {Number} beta The mean packing ratio [-]
   *
   * @return {Number} xi
   */
  static propagatingFluxRatio(sigma, beta) {
    if (sigma <= 0) {
      throw new RangeError("sigma must be > 0.")
    }
    const xi = Math.exp((0.792 + 0.681 * Math.sqrt(sigma)) * (beta + 0.1)) / (192 + 0.2595 * sigma)
    return xi
  }

  /**
   * Calculates the effective heating number: epsilon.
   *
   * Rothermel 1972: eq. (14) and (77).
   *
   * @param {Number} sv The SAV ratio value for an individual particle [ft2/ft3].
   *
   * @return {Number} epsilon.
   */
  static effectiveHeatingNumber(sv) {
    let epsilon = 0
    if (sv > 0) {
      epsilon = Math.exp(-138 / sv)
    }
    return epsilon
  }

  /**
   * Calculates the heat of preignition: Q_ig.
   *
   * Rothermel 1972: eq. (12) and (78).
   *
   * @param {Number} Mf The fuel moisture value for an individual fuel particle [%].
   *
   * @return {Number} Q_ig.
   */
  static heatOfPreignition(Mf) {
    const Q_ig = 250 + 1116 * (Mf * 0.01) // Mf = [fraction]
    return Q_ig
  }

  /**
   * Calculates the heat sink term: hsk.
   *
   * Rothermel 1972: eq. (77).
   *
   * @param {Number[]} preignitionHeat An array of heat of preignition values for individual particles
   * (Q_ig).
   * @param {Number[]} effectiveHeating An array of effective heating number values for individual particles
   * (epsilon).
   * @param {Number[]} sw An array of (sv * w0) weighting values for individual fuel particles (sw).
   * @param {Number} density The mean bulk density for the fuel complex (rho_b).
   *
   * @return {NUmber} hsk [Btu/ft3]
   */
  static heatSink(preignitionHeat, effectiveHeating, sw, density) {
    let Qig_t = 0   // sum[i=1,n][Qig_i]
    let sw_t = 0    // sum[i=1,n][sw_i]
    const len = sw.length
    for (let i = 0; i < sw.length; i++) {
      Qig_t += preignitionHeat[i] * effectiveHeating[i] * sw[i];
      sw_t += sw[i];
    }
    const hsk = density * (Qig_t / sw_t);
    return hsk;
  }

  /**
   * Calculates the wind factor: phi_w.
   *
   * Rothermel 1972: eq. (47) and (79),(82),(83),(84)
   *
   * @param {Number} midFlameWindSpd The wind speed at mid-flame height [ft/min].
   * @param {Number} sigma The characteristic SAV ratio for the fuelbed [ft2/ft3].
   * @param {Number} beta_ratio The relative packing ratio [beta/beta_opt].
   *
   * @return {Number} phi_w
   */
  static windFactor(midFlameWindSpd, sigma, beta_ratio) {
    const C = Rothermel.windParameterC(sigma)
    const B = Rothermel.windParameterB(sigma)
    const E = Rothermel.windParameterE(sigma)
    const phi_w = Rothermel.windFactor2(midFlameWindSpd, C, B, E, beta_ratio)
    return phi_w
  }

  /**
   * Calculates the wind multiplier for the rate of spread: phi_w.
   *
   * Rothermel 1972: eq. (47)
   *
   * @param {Number} midFlameWindSpd The wind speed at mid-flame height [ft/min].
   * @param {Number} C Result from Rothermel 1972: eq. (48).
   * @param {Number} B Result from Rothermel 1972: eq. (49).
   * @param {Number} E Result from Rothermel 1972: eq. (50).
   * @param {Number} beta_ratio The relative packing ratio [beta/beta_opt].
   *
   * @return {Number} phi_w
   */
  static windFactor2(midFlameWindSpd, C, B, E, beta_ratio) {
    const phi_w = C * Math.pow(midFlameWindSpd, B) * Math.pow(beta_ratio, -E)
    return phi_w
  }

  /**
   * Calculates the wind parameter C.
   *
   * Rothermel 1972: eq. (48)
   *
   * @param {Number} sigma The characteristic SAV ratio for the fuelbed [ft2/ft3].
   *
   * @return {Number} C
   */
  static windParameterC(sigma) {
    const C = 7.47 * Math.exp(-0.133 * Math.pow(sigma, 0.55))
    return C
  }

  /**
   * Calculates the wind parameter B.
   *
   * Rothermel 1972: eq. (49)
   *
   * @param {Number} sigma The characteristic SAV ratio for the fuelbed [ft2/ft3].
   *
   * @return {Number} B
   */
  static windParameterB(sigma) {
    const B = 0.02526 * Math.pow(sigma, 0.54)
    return B
  }

  /**
   * Calculates the wind parameter E.
   *
   * Rothermel 1972: eq. (50)
   *
   * @param {Number} sigma The characteristic SAV ratio for the fuelbed [ft2/ft3].
   *
   * @return {Number} E
   */
  static windParameterE(sigma) {
    const E = 0.715 * Math.exp(-0.000359 * sigma)
    return E;
  }

  /**
   * Calculates the slope multiplier for the rate of spread: phi_s.
   *
   * Rothermel 1972: eq. (51) and (78)
   *
   * @param {Number} slopeDegrees The steepness of the slope [degrees].
   * @param {Number}bbeta The mean packing ratio.
   *
   * @return {Number} phi_s
   */
  static slopeFactor(slopeDegrees, beta) {
    const phi = slopeDegrees * TO_RADIANS
    const tan_phi = Math.tan(phi)
    const phi_s = 5.275 * Math.pow(beta, -0.3) * Math.pow(tan_phi, 2)
    return phi_s
  }

  /**
   * Calculates the effective wind speed from the combined wind and slope factors: efw.
   *
   * Rothermel 1972: eq. (87)
   *
   * @param {NUmber} phiEw The combined wind and slope factors [phiW + phiS].
   * @param {NUmber} beta_ratio beta/beta_opt.
   * @param {NUmber} sigma The characteristic SAV ratio [ft2/ft3].
   * 
   * @return {NUmber} efw [ft/min]
   */
  static effectiveWindSpeed(phiEw, beta_ratio, sigma) {
    const C = Rothermel.windParameterC(sigma)
    const B = Rothermel.windParameterB(sigma)
    const E = Rothermel.windParameterE(sigma)
    const efw = Rothermel.effectiveWindSpeed2(phiEw, C, B, E, beta_ratio)
    return efw
  }

  /**
   * Calculates the effective wind speed from the combined wind and slope factors: efw.
   *
   * Rothermel 1972: eq. (87)
   *
   * @param {NUmber} phiEw The combined wind and slope factors [phiW + phiS].
   * @param {NUmber} C Result from Rothermel 1972: eq. (48).
   * @param {NUmber} B Result from Rothermel 1972: eq. (49).
   * @param {NUmber} E Result from Rothermel 1972: eq. (50).
   * @param {NUmber} beta_ratio
   * 
   * @return {NUmber} efw [ft/min]
   */
  static effectiveWindSpeed2(phiEw, C, B, E, beta_ratio) {
    // Effective windspeed: actually this is only the inverse function of phi_w
    const efw = (Math.pow(phiEw / (C * Math.pow(beta_ratio, -E)), 1 / B))
    return efw
  }

  /**
   * Calculates the rate of spread with wind and/or slope: ros.
   *
   * Rothermel 1972: eq. (52) - heat source / heat sink
   *
   * @param {Number} reactionIntensity The fire reaction intensity (I_r) [BTU/ft2/min].
   * @param {Number} propogatingFlux The fire propagating flux (xi) [fraction].
   * @param {Number} windFactor The wind coefficient (phi_w).
   * @param {Number} slopeFactor The slope coefficient (phi_s).
   * @param {Number} heatSink The total heat sink (hsk) [Btu/ft3].
   *
   * @return {Number} ros [ft/min]
   */
  static rateOfSpread(reactionIntensity, propogatingFlux, windFactor, slopeFactor, heatSink) {
    if (heatSink <= 0) {
      throw new RangeError("heatSink must be > 0.");
    }
    const ros = (reactionIntensity * propogatingFlux * (1 + windFactor + slopeFactor)) / heatSink;
    return ros;
  }

  /**
   * Calculates the rate of spread without wind and slope: ros.
   *
   * Rothermel 1972: eq. (52) - heat source / heat sink
   *
   * @param {Number} reactionIntensity The fire reaction intensity (I_r) [BTU/ft2/min].
   * @param {Number} propogatingFlux The fire propagating flux (xi) [fraction].
   * @param {Number} heatSink The total heat sink (hsk) [Btu/ft3].
   *
   * @return {Number} ros [ft/min]
   */
  static rateOfSpreadNoWindNoSlope(reactionIntensity, propogatingFlux, heatSink) {
    if (heatSink <= 0) {
      throw new RangeError("heatSink must be > 0.")
    }
    const ros = (reactionIntensity * propogatingFlux) / heatSink
    return ros
  }

  /**
   * Calculates the flame zone depth: fzd.
   *
   * The depth, or front-to-back distance, of the actively flaming zone of a free spreading fire
   * can be determined from the rate of spread and the particle-residence time. Albini 1976: pg.
   * 86
   *
   * @param {Number} rateOfSpread The fire rate of spread (ros) [ft/min].
   * @param {Number} flameResidenceTime The fuelbed's flame residence time [min].
   *
   * @return {Number} fzd [ft]
   */
  static flameZoneDepth(rateOfSpread, flameResidenceTime) {
    const fzd = rateOfSpread * flameResidenceTime
    return fzd
  }

  /**
   * Calculates Byram's fireline intensity: I.
   *
   * Byram's intensity, I, is the rate of heat release per unit of fire edge. The reaction
   * intensity, I_r, provided by Rothermel's spread model is the rate of energy release per unit
   * area in the actively flaming zone.
   *
   * Albini 1976: eq. (16), pg. 86
   *
   * @param flameZoneDepth The depth of the actively flaming zone [ft].
   * @param reactionIntensity The fuelbed's fire reaction intensity (I_r) [Btu/ft2/min].
   *
   * @return I [Btu/ft/s]
   */
  static firelineIntensity(flameZoneDepth, reactionIntensity) {
    const I = reactionIntensity * flameZoneDepth / 60
    return I
  }

  /**
   * Calculates flame length: L.
   *
   * Albini 1976: eq. (17) pg. 86
   *
   * @param {Number} firelineIntensity Byram's fireline intensity (I) [Btu/ft/s].
   *
   * @return {Number} L [ft]
   */
  static flameLength(firelineIntensity) {
    const L = 0.45 * Math.pow(firelineIntensity, 0.46)
    return L;
  }

  /**
   * Calculates the wind adjustment factor for scaling wind speed from 20-ft to midflame height.
   *
   * Wind adjustment factor is calculated as an average from the top of the fuel bed to twice the
   * fuel bed depth, using Albini and Baughman (1979) equation 9 (page 5).
   *
   * @param {Number} fuelDepth Fuel bed depth (height) [ft].
   * 
   * @return {Number} Wind adjustment factor, waf [0..1]
   */
  static midFlameWindAdjustmentFactor(fuelDepth) {
    let waf = 1.0
    if (fuelDepth > 0) {
      // From BehavePlus5, xfblib.cpp by Collin D. Bevins
      waf = 1.83 / Math.log((20 + 0.36 * fuelDepth) / (0.13 * fuelDepth))
    }
    return Math.min(Math.max(waf, 0), 1)
  }

  /**
   * Computes the mid-flame wind speed from 20 foot wind speeds. Used to compute rate of spread.
   *
   * @param {Number} wndSpd20Ft Wind speed 20 feet above the vegetation [MPH]
   * @param {Number} fuelDepth Vegetation height [feet]
   * 
   * @return {Number} Wind speed at vegetation height
   */
  static calcWindSpeedMidFlame(wndSpd20Ft, fuelDepth) {
    return wndSpd20Ft * Rothermel.midFlameWindAdjustmentFactor(fuelDepth)
  }

  /**
   * Computes the wind speed at the fuel level from 20 foot wind speeds. Used to compute wind
   * cooling effect on fuel temperatures.
   *
   * @param {Number} wndSpd20Ft Wind speed 20 feet above the vegetation [MPH]
   * @param {Number} fuelDepth Vegetation height [feet]
   * 
   * @return {Number} Wind speed at vegetation height
   */
  static calcWindSpeedNearFuel(wndSpd20Ft, fuelDepth) {
    // Equation #36
    // The ratio of windspeed at vegetation height to that at
    // 20 feet above the vegitation is given by:
    //  U_h' / U_20+h' = 1 / ln((20 + 0.36 * h') / 0.13 * h')
    //      where:
    //          h' = vegitation height [feet]
    if (fuelDepth == 0) {
      fuelDepth = 0.1
    }
    const U_h = (1.0 / Math.log((20 + 0.36 * fuelDepth) / (0.13 * fuelDepth))) * wndSpd20Ft
    return U_h
  }

  /**
   * Calculates the fire ellipse eccentricity from the effective wind speed.
   *
   * Anderson 1983: eq. (4)
   *
   * <pre>
   * Consider using Anderson 1983: eq. (17)
   *      l/w = 0.936 EXP(0.1l47U) + 0.461 EXP(-0.0692U)
   *  where U = midflame miles per hour.
   * </pre>
   * 
   * @param {Number} effectiveWind The effective wind speed of the combined wind and slope. [mph]
   * 
   * @return {Number} The eccentricity of the ellipse
   */
  static eccentricity(effectiveWind) {
    let eccentricity = 0
    if (effectiveWind > 0) {
      // From FireLib 1.04, firelib.c by Collin D. Bevins
      // a1 = major axis of semiellipse at the rear of the fire,
      //    = 1. + 0.25 * effectiveWindSpd / 88.0);
      const lbRatio = 1. + 0.002840909 * effectiveWind
      if (lbRatio > 1.00001) {
        eccentricity = Math.sqrt(Math.pow(lbRatio, 2) - 1.0) / lbRatio
      }
    }
    return eccentricity
  }

  /**
   * Compute difference between fuel temperature and the air temperature due to solar heating and
   * wind cooling effects.
   *
   * Rothermel 1986: eq. (1)
   *
   * @param {Number} I radiation intensity [cal/c2 * min]
   * @param {Number} T_a temperature of air [fahrenheit]
   * @param {Number} U_h wind velocity at fuel level [mph]
   *
   * @return {Number} T_f temperature of fuel [fahrenheit]
   */
  static calcFuelTemp(I, T_a, U_h) {
    // Rothermel et al, 1986, page 9
    // Equation #1
    // The difference in temperature between the air and fuel is assumed
    // to be directly proportional to the incident radiation intensity, I,
    // and inversely proportional to the wind velocity, U, and two constants
    // attributed to fuel conditions:
    //  T_f - T_a == I / (0.015 * U_h + 0.026)
    // where:
    //  T_f  =   temperature of fuel [farenheit]
    //  T_a  =   temperature of air  [farenheit]
    //  I    =   radiation intencity [cal/c2 * min]
    //  U_h  =   wind velocity at fuel level [mph]

    const T_f = (I / (0.015 * U_h + 0.026)) + T_a;
    return T_f;
  }

  /**
   * Compute the relative humidity for the air immediately adjacent to the fuel.
   *
   * Rothermel 1986: eq. (2)
   *
   * @param {Nunber} H_a relative humidity of the air [percent]
   * @param {Nunber} T_f fuel temperature [Fahrenheit]
   * @param {Nunber} T_a air temperature [Fahrenheit]
   * 
   * @return {Nunber} H_f - relative humidity of the air next to the fuel [percent]
   */
  static calcRelativeHumidityNearFuel(H_a, T_f, T_a) {
    // Rothermel et al, 1986, page 9
    // Equation #2
    // Correction for relative humidity as a function of the fuel temperature
    // and air temperature:
    //  H_f = H_a * exp(-0.033(T_f - T_a))
    const H_f = H_a * Math.exp(-0.033 * (T_f - T_a))
    return H_f
  }

  /**
   * Computes the earth-sun (center of mass) distance squared.
   *
   * @param {Number} delta solar declination [radians]
   * 
   * @return {Nunber} r2 earth-sun distance squared
   */
  static calcEarthSunDistanceSqrd(delta) {
    // Rothermel et al, 1986, page 11
    // Equation #7
    // The earth-sun distance (squared) by analytic solution to tabular values
    const r2 = 0.999847 + (0.001406 * (delta * TO_DEGREES));
    return r2;
  }

  /**
   * Computes the solar irradiance.
   *
   * Rothermel 1986: eq. (3)
   *
   * @param {Nunber} I_a irradiance at the forest floor perpendicular to the solar ray [cal/cm2*min]
   * @param {Nunber} r2 The earth-sun (center of mass) distance squared
   * @param {Nunber} A solar elevation angle to the sun [radians]
   *
   * @return {Nunber} I - incident radiation on the forest floor [cal/cm2*min]
   */
  static calcSolarIrradianceOnHorzSurface(I_a, r2, A) {
    // Rothermel et al, 1986, page 9
    // Equation #3
    // I = (I_a / r2) * sin A
    if (A <= 0) {
      return 0
    }
    const I = (I_a / r2) * Math.sin(A);
    return I;
  }

  /**
   * Computes the optical air mass, i.e., the ratio of the optical path length of radiation
   * through the atmosphere at angle A, to the path length toward the zenith at sea level.
   *
   * Rothermel 1986: eq. (16)
   *
   * @param {Number} A the solar altitude angle [radians]
   * @param {Number} E the elevation at angle A [feet]
   * 
   * @return {Number} M the optical air mass ratio
   */
  static calcOpticalAirMass(A, E) {
    // Equation #16
    // M = (absolute_pressure / sea_level_pressure) 
    // csc A = exp(-0.0000448E) csc A
    let M = 0;
    if (A > 0) {
      M = Math.exp(-0.0000448 * E) * 1.0 / Math.sin(A);
    }
    return M;
  }

  /**
   * Computes irradiance at forest floor perpendicular to solar ray (1 [cal/cm2*min] = 697.8
   * [watts/m2])
   *
   * @param {Number} M is the optical air mass ratio
   * @param {Number} S_c cloud cover [percent]
   * @param {Number} p is the transparency coefficient
   * 
   * @return {Number} attenuated irradiance [cal/cm2*min]
   */
  static calcAttenuatedIrradiance(M, S_c, p) {
    // I_a = I_M * tau_n
    //  where:
    //      tau_n is net transmittance of clouds and trees
    //      I_M is the direct solar irradiance including atmospheric attenuation
    // I_M = I_o * pM
    //  where:
    //      I_o is incident intensity or solar constant 1.98 cal/cm2*min
    //      p is the transparency coefficient
    //      M is the optical air mass, the ratio of the
    //
    if (M <= 0) {
      return 0;
    }
    const I_o = 1.98;                 // solar constant [cal/cm2*min]
    const tau_c = 1 - (S_c / 100.0);  // cloud shade transmittance
    const tau_t = 1;                  // tree shade transmittance (default 1 for now)
    const tau_n = tau_t * tau_c;
    const I_M = I_o * Math.pow(p, M); // incident radiation attenuated by atmosphere
    const I_a = I_M * tau_n;          // irradiance at forest floor perpendicular to solar ray
    return I_a;
  }

  /**
   * Computes the irradiance on a slope (neglecting the small variation in r).
   *
   * Rothermel 1986, eq. (9),(10) and (11), adjusted for Z relative to North instead of East
   *
   * @param {Number} alpha slope angle from horizontal at slope azimuth [radians]
   * @param {Number} beta aspect of the slope [radians]
   * @param {Number} A solar altitude (elevation) angle [radians]
   * @param {Number} Z solar azimuth angle (true) [radians]
   * @param {Number} I_a attenuated irradiance [cal/cm2*min]
   * 
   * @return {Number} incident radiation intensity [cal/cm2*min]
   */
  static calcIrradianceOnASlope(alpha, beta, A, Z, I_a) {
    // Rothermel et al, 1986, page 11
    // Equation #9, 10 and 11
    //  I = Ia * sin zeta
    // where:
    //  alpha   = slope angle from horizontal at slope azimuth
    //  psi     = slope angle at solar azimuth Z
    //      zeta replaces A in equation #3, the solar angle to the slope in
    //      the plane normal to the slope
    //  sin zeta = sin(A - psi) * cos(alpha) / cos(psi)
    //  tan psi  = tan alpha * sin(z - beta)
    // where:
    //      (A - psi) is the solar angle to the slope in local vertical plane
    //      and psi is the slope angle at the solar azimuth, z,
    //

    // Precondition: Sun above the horizon
    if (A <= 0) {
      return 0
    }
    // Precondition: Must have sunlight (not total shade)
    if (I_a <= 0) {
      return 0
    }

    // Adjusted original algorithm from East azimuth to North azimuth with 1/2 PI.
    const tanPsi = Math.tan(alpha) * Math.sin(Z - beta - HALF_PI)
    const psi = Math.atan(tanPsi)
    const sinZeta = Math.sin(A - psi) * Math.cos(alpha) / Math.cos(psi)

    // I is the incident radiation intensity
    const I = I_a * sinZeta

    // Post condition: I >= 0
    return (I > 0) ? I : 0
  }

  /**
   * Computes the Canadian Standard Daily Fine Fuel Moisture Code (FFMC) and converts it to a
   * percentage. Calculates the fuel moisture percentage for early afternoon from noon-time
   * weather conditions or forecasts. <br/>
   *
   * Rothermel, et al, "Modeling moisture content of fine dead wildland fuels: input to the BEHAVE
   * fire prediction system." Research Paper INT-359. 1986. Equations located on page 47. <br/>
   *
   * Note: FFMC: Low = 0 - 72; Moderate = 73 - 77; High = 78 - 82; Extreme > 82
   *
   * @param {Number} m_0 initial fuel moisture at noon [percent]
   * @param {Number} T_f air temp immediately adjacent to fuel [Fahrenheit]
   * @param {Number} H_f relative humidity immediately adjacent to fuel [percent]
   * @param {Number} W 20 foot wind speed [MPH]
   * @param {Number} R rainfall amount [inches]
   * 
   * @return {Number} Fuel moisture percent derived from computed FFMC code [percent]
   */
  static calcCanadianStandardDailyFineFuelMoisture(m_0, T_f, H_f, W, R) {

    // f_0 - initial moisture converted to a FFMC
    const f_0 = 101 - m_0

    // Equation #1 - adust the initial fuel moisture code (f0) for rain
    let f_R = 0     // f_0 modified for rain
    // if rain > 0.02"
    if (R > 0.02) {
      const R_A = Math.min(R, 1.5)
      let F = 0
      if (R_A <= 0.055) { // [inches]
        F = -56.0 - 55.6 * Math.log(R_A + 0.04)
      } else if (R_A <= 0.225) { // [inches]
        F = -1.0 - 18.2 * Math.log(R_A - 0.04)
      } else {
        F = 14 - 8.25 * Math.log(R_A - 0.075)
      }
      f_R = Math.max(0, (F * f_0 / 100) + 1 - 8.73 * Math.exp(-0.1117 * f_0))
    } else {
      // little or no rain
      f_R = f_0
    }

    // Equation #2
    // m_R - initial fuel moisture [percent] adjusted for rain
    const m_R = 101 - f_R

    // Equation #3 - equilibrium drying curve
    const E_d = (0.942 * Math.pow(H_f, 0.679)) + 11 * Math.exp((H_f / 10) - 10);

    // Equation #4 - equilibrium wetting curve
    const E_w = (0.597 * Math.pow(H_f, 0.768)) + 14 * Math.exp((H_f / 8) - 12.5);

    // m - fine fuel moisture adjusted for humidity and wind
    let m = 0;
    if (Math.abs(m_R - E_d) < 0.0001) { // nearly equal
      m = m_R;
    } else if (m_R < E_d) {
      // Wetting: fuel moisture is below the drying curve so a wetting trend is in effect
      // Equation #5
      m = E_w + (m_R - E_w) / 1.9953;       //-- original
      //m = E_W + (E_W - m_R) / 1.9953;     // corrected based on Anderson 2009 87-10 
    } else { 
      // Drying: fuel moisture is above the drying curve so a drying trend is in effect
      // Here we constrain 20' wind to between 1 and 14 mph
      W = Math.min(Math.max(1.0, W), 14.0);
      // Equations #6 and #7
      const x = 0.424 * (1 - Math.pow(H_f / 100, 1.7)) + 0.088 * Math.pow(W, 0.5) * (1 - Math.pow(H_f / 100, 8));
      m = E_d + (m_R - E_d) / Math.pow(10, x);
    }
    // compute fine fuel moisture delta for temperature
    let delta = 0;
    if (f_0 < 99.0) {
      delta = Math.max(-16.0, (T_f - 70) * (0.63 - 0.0065 * f_R));
    }
    // final FFMC code constrained to between 0 and 99
    const f = Math.max(0, Math.min(99, 101 - m + delta));

    // FFMC code converted to fuel moisture percentage
    return 101 - f;
  }

//    /**
//     * Computes the Canadian Hourly Fine Fuel Moisture percentage.
//     *
//     * Note: The calcCanadianStandardDailyFineFuelMoisure provides the first m_0. Subsequently, the
//     * previous hour's m becomes m_0.
//     *
//     * Rothermel, et al, "Modeling moisture content of fine dead wildland fuels: input to the BEHAVE
//     * fire prediction system." Research Paper INT-359. 1986. Equations located on page 47.
//     *
//     * @param m_0 previous hour's fuel moisture [percent]
//     * @param H relative humidity [percent]
//     * @param T_c air temperature [Celsius]
//     * @param W_k 20' wind speed [KPH]
//     *
//     * @return fuel moisture [percent]
//     */
//    static public double calcCanadianHourlyFineFuelMoisture(double m_0, double H, double T_c, double W_k) {
//        // Equation #1 [not used/applicatble] converts previous hours FFMC to fuel moisture percent
//        // m_0 = previous hour's fuel moisture;
//
//        // constrain wind to 22.5 kph
//        W_k = min(max(0, W_k), 22.5);
//
//        double m = 0;
//
//        // Equation #2a compute equilibruim moisture curve (EMC) for drying
//        double E_d = 0.942 * pow(H, 0.679) + 11 * exp((H - 100) / 10)
//                + 0.18 * (21.1 - T_c) * (1 - exp(-0.115 * H));
//        // Equation #2b compute equilibruim moisture curve (EMC) for wetting
//        double E_w = 0.618 * pow(H, 0.753) + 10 * exp((H - 100) / 10)
//                + 0.18 * (21.1 - T_c) * (1 - exp(-0.115 * H));
//
//        if (m_0 > E_d) {
//            // Drying...
//            // Equations #3a and #3b compute log drying rate for hourly computation, log base 10
//            double k_a = 0.424 * (1 - pow(H / 100, 1.7)) + 0.0694 * pow(W_k, 0.5) * (1 - pow(H / 100, 8));
//            double k_d = k_a * 0.0579 * exp(0.0365 * T_c);
//            // Equation #5a computes final fuel moisture percent
//            m = E_d + (m_0 - E_d) * exp(-2.303 * k_d);
//
//        } else if (m_0 < E_w) {
//            // Wetting...
//
//            // Rothermel (4a) and (4b) compute log wetting rate for hourly computation, log base 10
//            double k_b = 0.424 * (1 - (pow((100 - H) / 100, 1.7)))
//                    + 0.0694 * pow(W_k, 0.5) * (1 - (pow((100 - H) / 100, 8)));
//
//            double k_w = k_b * 0.0579 * exp(0.0365 * T_c);
//
//            // Equation #5b computes final fuel moisture percent
//            m = E_w - (E_w - m_0) * exp(-2.303 * k_w);    // Rothemel pg 48
//
//        } else {
//            m = m_0;
//        }
//        return m;
//    }
//
//    /**
//     * Computes the fine dead fuel moisture using the EMC from Canadian Standard Daily Fine Fuel
//     * Moisture Code formula. The return value assumes instantaneous drying and wetting. <br/>
//     *
//     * Anderson, K., 2009, "A Comparison of Hourly Fire Fuel Moisture Code Calculations within
//     * Canada", Canadian Forest Service
//     *
//     * @param m_0 Initial fuel moisture used to determine drying or wetting phase. [percent]
//     * @param T_c Air temperature immediately adjacent to fuel [Celsius]
//     * @param H Relative humidity immediately adjacent to fuel [percent]
//     * @return Dead fine fuel moisture [percent]
//     */
//    static public double calcFineDeadFuelMoisture(double m_0, double T_c, double H) {
//
//        // Van Wagner Eq. #2a (87-8a) equilibruim moisture curve (EMC) for drying
//        double E_d = 0.942 * pow(H, 0.679) + 11 * exp((H - 100) / 10)
//                + 0.18 * (21.1 - T_c) * (1 - exp(-0.115 * H));
//
//        // Van Wagner Eq. #2b (87-8b) equilibruim moisture curve (EMC) for wetting
//        double E_w = 0.618 * pow(H, 0.753) + 10 * exp((H - 100) / 10)
//                + 0.18 * (21.1 - T_c) * (1 - exp(-0.115 * H));
//
//        // m - fine fuel moisture 
//        double m;
//        if (m_0 > E_d) {
//            // Instantaneous Drying
//            m = E_d;
//        } else if (m_0 < E_w) {
//            // Instantaneous Wetting
//            m = E_w;
//        } else {
//            // No change
//            m = m_0;
//        }
//        return max(m, 0);
//    }

}