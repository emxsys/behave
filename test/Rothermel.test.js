
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