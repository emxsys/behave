
import Rothermel from '../src/Rothermel'

describe('Rothermel', () => {
  it('should exist', () => {
    expect(Rothermel).toBeTruthy()
  })
})


describe('meanBulkDensity', () => {
  const w0 = [1, 2, 3]  // sum = 6
  const height = 2
  const rho_b = 3

  // expect().toThrowError requires the function to be tested be wrapped in an anonymous function
  it('throws with height zero', () => {
    expect(() => {
      Rothermel.meanBulkDensity(w0, 0)
    }).toThrowError(RangeError)
  })
  
  it('throws with height < zero', () => {
    expect(() => {
      Rothermel.meanBulkDensity(w0, -1)
    }).toThrowError(RangeError)
  })

  it('result is valid', () => {
    expect(Rothermel.meanBulkDensity(w0, height)).toBe(rho_b)
  })
})


describe('meanPackingRatio', () => {
  const rho_b = 1
  const rho_p = 10
  const beta = 0.1

  // Note: expect().toThrowError requires the function to be tested be wrapped in an anonymous function
  
  it('throws with rho_p zero', () => {
    expect(() => {
      Rothermel.meanPackingRatio(rho_b, 0)
    }).toThrowError(RangeError)
  })
  
  it('throws with rho_p < zero', () => {
    expect(() => {
      Rothermel.meanPackingRatio(rho_b, -1)
    }).toThrowError(RangeError)
  })
  
  it('throws with beta < zero', () => {
    expect(() => {
      Rothermel.meanPackingRatio(rho_b, -1)
    }).toThrowError(RangeError)
  })
  
  it('throws with beta > 0.12', () => {
    expect(() => {
      Rothermel.meanPackingRatio(rho_b, 7.6)
    }).toThrowError(RangeError)
  })
  
  it('throws with beta < 0', () => {
    expect(() => {
      Rothermel.meanPackingRatio(-1, rho_p)
    }).toThrowError(RangeError)
  })

  it('result is valid', () => {
    expect(Rothermel.meanPackingRatio(rho_b, rho_p)).toBe(beta)
  })

})


describe('optimalPackingRatio', () => {
    const sigma = 1
    const beta_opt = 3.348 * Math.pow(sigma, -0.8189)

  it('algorithm has not changed', () => {
    expect(Rothermel.optimalPackingRatio(sigma)).toBe(beta_opt)
  })
})


describe('characteristicSAV', () => {
    const sv = [3]  // ft2/ft3
    const w0 = [2]  // lbs/ft2
    const sigma = (sv[0] * sv[0] * w0[0]) / (sv[0] * w0[0])

  it('algorithm has not changed', () => {
    expect(Rothermel.characteristicSAV(sv, w0)).toBe(sigma)
  })
})


describe('reactionVelocity', () => {
    const sigma = 2 
    const beta_ratio = 3
    const gamma = 4.702292202212657e-33

  it('algorithm has not changed', () => {
    expect(Rothermel.reactionVelocity(sigma, beta_ratio)).toBe(gamma)
  })
})


describe('reactionIntensity', () => {
    const gamma = 2 
    const heat = 3
    const eta_M = 4
    const eta_s = 5
    const I_r = gamma * heat * eta_M * eta_s

  it('algorithm has not changed', () => {
    expect(Rothermel.reactionIntensity(gamma, heat, eta_M, eta_s)).toBe(I_r)
  })
})


describe('flameResidenceTime', () => {
    const sigma = 2 
    const tau = 384 / sigma

  it('algorithm has not changed', () => {
    expect(Rothermel.flameResidenceTime(sigma)).toBe(tau)
  })
})


describe('heatRelease', () => {
    const I_r = 2 
    const tau = 3
    const hpa = I_r * tau

  it('algorithm has not changed', () => {
    expect(Rothermel.heatRelease(I_r, tau)).toBe(hpa)
  })
})