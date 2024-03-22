#ifndef G_RANDOM_ENGINE_H
#define G_RANDOM_ENGINE_H

namespace Garfield {

/// Abstract base class for random number generators.

class RandomEngine {
 public:
  /// Constructor
  RandomEngine() = default;
  /// Destructor
  virtual ~RandomEngine() {}

  /// Draw a random number.
  virtual double Draw() = 0;
  /// Initialise the random number generator.
  virtual void Seed(const unsigned int s) = 0;
  /// Print some information about the random number generator.
  virtual void Print() = 0;
};
}

#endif
