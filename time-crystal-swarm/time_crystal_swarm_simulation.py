"""
Time Crystal Swarm Simulation

Simulates distributed energy extraction from coupled time crystal systems.
Demonstrates swarm stability under weak measurement.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, List
import warnings
warnings.filterwarnings('ignore')


class TimeCrystalSwarm:
    """Simulate a swarm of coupled time crystals with weak measurement"""
    
    def __init__(self, n_crystals: int = 25, coupling_strength: float = 0.1):
        self.n_crystals = n_crystals
        self.coupling = coupling_strength
        self.n_core = int(n_crystals * 0.4)
        
        self.states = None
        self.coherence_history = []
        self.energy_history = []
        self.extracted_energy = []
        
        print(f"Initialized swarm: {n_crystals} crystals")
        print(f"Core: {self.n_core}, Peripheral: {n_crystals - self.n_core}")
        
    def initialize_states(self):
        """Initialize all crystals in superposition state"""
        self.states = np.ones((self.n_crystals, 2)) / np.sqrt(2)
        
    def time_crystal_evolution(self, crystal_idx: int, time: float) -> np.ndarray:
        """Evolve time crystal state with period doubling"""
        omega = 2 * np.pi
        phase = omega * time / 2
        
        cos_half = np.cos(phase / 2)
        sin_half = np.sin(phase / 2)
        
        rotation = np.array([
            [cos_half, -sin_half],
            [sin_half, cos_half]
        ])
        
        return rotation @ self.states[crystal_idx]
    
    def apply_coupling(self):
        """Apply inter-crystal coupling to create entanglement"""
        new_states = self.states.copy()
        
        for i in range(self.n_crystals):
            left = (i - 1) % self.n_crystals
            right = (i + 1) % self.n_crystals
            
            neighbor_avg = (self.states[left] + self.states[right]) / 2
            new_states[i] = (1 - self.coupling) * self.states[i] + \
                           self.coupling * neighbor_avg
            
            norm = np.linalg.norm(new_states[i])
            if norm > 0:
                new_states[i] /= norm
        
        self.states = new_states
    
    def weak_measurement(self, crystal_idx: int, strength: float = 0.01) -> float:
        """Perform weak measurement on peripheral crystal"""
        prob_0 = np.abs(self.states[crystal_idx, 0])**2
        prob_1 = np.abs(self.states[crystal_idx, 1])**2
        
        outcome = np.random.rand() < prob_0
        
        if outcome:
            self.states[crystal_idx, 0] += strength
            self.states[crystal_idx, 1] -= strength * 0.5
        else:
            self.states[crystal_idx, 1] += strength
            self.states[crystal_idx, 0] -= strength * 0.5
        
        norm = np.linalg.norm(self.states[crystal_idx])
        self.states[crystal_idx] /= norm
        
        energy = strength * (prob_1 - prob_0)
        return energy
    
    def calculate_coherence(self, crystal_indices: List[int]) -> float:
        """Calculate average coherence of specified crystals"""
        coherence_sum = 0
        
        for idx in crystal_indices:
            purity = np.abs(np.sum(self.states[idx]**2))
            coherence_sum += purity
        
        return coherence_sum / len(crystal_indices)
    
    def calculate_total_energy(self) -> float:
        """Calculate total energy in swarm"""
        energy = 0
        for i in range(self.n_crystals):
            energy += np.abs(self.states[i, 1])**2
        return energy
    
    def simulate(self, n_steps: int = 1000, measurement_strength: float = 0.01):
        """Run swarm simulation"""
        print(f"\nRunning simulation: {n_steps} steps, epsilon={measurement_strength}")
        
        self.initialize_states()
        
        core_indices = list(range(self.n_core))
        peripheral_indices = list(range(self.n_core, self.n_crystals))
        
        for step in range(n_steps):
            for i in range(self.n_crystals):
                self.states[i] = self.time_crystal_evolution(i, step * 0.1)
            
            self.apply_coupling()
            
            total_extracted = 0
            for i in peripheral_indices:
                if np.random.rand() < 0.1:
                    energy = self.weak_measurement(i, measurement_strength)
                    total_extracted += energy
            
            core_coherence = self.calculate_coherence(core_indices)
            total_energy = self.calculate_total_energy()
            
            self.coherence_history.append(core_coherence)
            self.energy_history.append(total_energy)
            self.extracted_energy.append(total_extracted)
            
            if step % 200 == 0:
                print(f"Step {step}: Core coherence = {core_coherence:.3f}")
        
        print("Simulation complete")
        self.print_results()
    
    def print_results(self):
        """Print simulation results"""
        avg_coherence = np.mean(self.coherence_history[-100:])
        total_extracted = np.sum(self.extracted_energy)
        
        print(f"\nResults:")
        print(f"Average core coherence: {avg_coherence:.3f}")
        print(f"Total energy extracted: {total_extracted:.3f}")
        print(f"Extraction efficiency: {(total_extracted/len(self.extracted_energy)):.4f} per step")
        
        if avg_coherence > 0.8:
            print("Status: Core remains coherent")
        elif avg_coherence > 0.5:
            print("Status: Partial decoherence")
        else:
            print("Status: Core collapsed")
    
    def plot_results(self, filename: str = None):
        """Visualize simulation results"""
        fig, axes = plt.subplots(2, 1, figsize=(10, 8))
        
        steps = range(len(self.coherence_history))
        
        axes[0].plot(steps, self.coherence_history, color='blue', linewidth=1.5)
        axes[0].axhline(y=0.8, color='green', linestyle='--', label='Success threshold')
        axes[0].set_xlabel('Time Steps')
        axes[0].set_ylabel('Core Coherence')
        axes[0].set_title(f'Time Crystal Swarm Stability (N={self.n_crystals})')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        cumulative = np.cumsum(self.extracted_energy)
        axes[1].plot(steps, cumulative, color='orange', linewidth=1.5)
        axes[1].set_xlabel('Time Steps')
        axes[1].set_ylabel('Cumulative Energy Extracted')
        axes[1].set_title('Energy Extraction Over Time')
        axes[1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if filename:
            plt.savefig(filename, dpi=150, bbox_inches='tight')
            print(f"Plot saved: {filename}")
        else:
            plt.show()


def compare_swarm_sizes():
    """Compare different swarm sizes to demonstrate scaling"""
    print("\nSwarm size comparison")
    
    sizes = [10, 25, 50, 100]
    results = []
    
    for n in sizes:
        print(f"\nTesting N={n}")
        swarm = TimeCrystalSwarm(n_crystals=n, coupling_strength=0.1)
        swarm.simulate(n_steps=500, measurement_strength=0.01)
        
        avg_coherence = np.mean(swarm.coherence_history[-100:])
        total_extracted = np.sum(swarm.extracted_energy)
        efficiency = total_extracted / 500
        
        results.append({
            'n_crystals': n,
            'coherence': avg_coherence,
            'extracted': total_extracted,
            'efficiency': efficiency
        })
    
    print("\nComparison results:")
    print(f"{'N':<10} {'Coherence':<15} {'Total Energy':<15} {'Efficiency'}")
    for r in results:
        print(f"{r['n_crystals']:<10} {r['coherence']:<15.3f} {r['extracted']:<15.3f} {r['efficiency']:.5f}")
    
    return results


def critical_threshold_experiment():
    """Find critical measurement strength"""
    print("\nCritical threshold experiment")
    
    n_crystals = 50
    strengths = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1]
    results = []
    
    for epsilon in strengths:
        print(f"\nTesting epsilon={epsilon}")
        swarm = TimeCrystalSwarm(n_crystals=n_crystals)
        swarm.simulate(n_steps=500, measurement_strength=epsilon)
        
        avg_coherence = np.mean(swarm.coherence_history[-100:])
        results.append((epsilon, avg_coherence))
    
    print("\nThreshold results:")
    print(f"{'Measurement':<15} {'Coherence':<15} {'Status'}")
    for eps, coh in results:
        status = "Stable" if coh > 0.8 else "Collapsed"
        print(f"{eps:<15.3f} {coh:<15.3f} {status}")


def main():
    """Run complete simulation suite"""
    print("Time Crystal Swarm Energy Extraction Simulation\n")
    
    print("Experiment 1: Basic swarm demonstration")
    swarm = TimeCrystalSwarm(n_crystals=25, coupling_strength=0.1)
    swarm.simulate(n_steps=1000, measurement_strength=0.01)
    swarm.plot_results('images/swarm_basic.png')
    
    print("\nExperiment 2: Scaling analysis")
    compare_swarm_sizes()
    
    print("\nExperiment 3: Critical threshold")
    critical_threshold_experiment()
    
    print("\nAll experiments complete")


if __name__ == "__main__":
    main()