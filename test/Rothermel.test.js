
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