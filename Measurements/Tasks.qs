// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT license.

namespace Quantum.Kata.Measurements {
    
    open Microsoft.Quantum.Primitive;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Extensions.Convert;
    open Microsoft.Quantum.Extensions.Math;
    
    
    //////////////////////////////////////////////////////////////////
    // Welcome!
    //////////////////////////////////////////////////////////////////
    
    // "Measurements" quantum kata is a series of exercises designed
    // to get you familiar with programming in Q#.
    // It covers the following topics:
    //  - single-qubit measurements,
    //  - discriminating orthogonal and nonorthogonal states.
    
    // Each task is wrapped in one operation preceded by the description of the task.
    // Each task (except tasks in which you have to write a test) has a unit test associated with it,
    // which initially fails. Your goal is to fill in the blank (marked with // ... comment)
    // with some Q# code to make the failing test pass.
    
    // The tasks are given in approximate order of increasing difficulty; harder ones are marked with asterisks.
    
    //////////////////////////////////////////////////////////////////
    // Part I. Single-Qubit Measurements
    //////////////////////////////////////////////////////////////////
    
    // Task 1.1. |0⟩ or |1⟩ ?
    // Input: a qubit which is guaranteed to be in |0⟩ or |1⟩ state.
    // Output: true if qubit was in |1⟩ state, or false if it was in |0⟩ state.
    // The state of the qubit at the end of the operation does not matter.
    operation IsQubitOne (q : Qubit) : Bool {
        // The operation M will measure a qubit in the Z basis (|0⟩ and |1⟩ basis)
        // and return Zero if the observed state was |0⟩ or One if the state was |1⟩.
        // To answer the question, you need to perform the measurement and check whether the result
        // equals One - either directly or using library function IsResultOne.
        //
        // Replace the returned expression with (M(q) == One)
        // Then rebuild the project and rerun the tests - T101_IsQubitOne_Test should now pass!

        return (M(q) == One);
    }
    
    
    // Task 1.2. Set qubit to |0⟩ state
    // Input: a qubit in an arbitrary state.
    // Goal:  change the state of the qubit to |0⟩.
    operation InitializeQubit (q : Qubit) : Unit {
        if(M(q) == One){
            X(q);
        }
    }
    
    
    // Task 1.3. |+⟩ or |-⟩ ?
    // Input: a qubit which is guaranteed to be in |+⟩ or |-⟩ state
    //        (|+⟩ = (|0⟩ + |1⟩) / sqrt(2), |-⟩ = (|0⟩ - |1⟩) / sqrt(2)).
    // Output: true if qubit was in |+⟩ state, or false if it was in |-⟩ state.
    // The state of the qubit at the end of the operation does not matter.
    operation IsQubitPlus (q : Qubit) : Bool {
        H(q);
        return (M(q) == Zero);
    }
    
    
    // Task 1.4. |A⟩ or |B⟩ ?
    // Inputs:
    //      1) angle alpha, in radians, represented as Double
    //      2) a qubit which is guaranteed to be in |A⟩ or |B⟩ state
    //         |A⟩ =   cos(alpha) * |0⟩ + sin(alpha) * |1⟩,
    //         |B⟩ = - sin(alpha) * |0⟩ + cos(alpha) * |1⟩.
    // Output: true if qubit was in |A⟩ state, or false if it was in |B⟩ state.
    // The state of the qubit at the end of the operation does not matter.
    operation IsQubitA (alpha : Double, q : Qubit) : Bool {
        Ry(-2.0 * alpha, q);

        return (M(q) == Zero);
    }
    
    
    // Task 1.5. |00⟩ or |11⟩ ?
    // Input: two qubits (stored in an array) which are guaranteed to be in |00⟩ or |11⟩ state.
    // Output: 0 if qubits were in |00⟩ state,
    //         1 if they were in |11⟩ state.
    // The state of the qubits at the end of the operation does not matter.
    operation ZeroZeroOrOneOne (qs : Qubit[]) : Int {
        // ...
        if(M(qs[0]) == Zero){
            return 0;
        } else {
            return 1;
        }
    }
    
    
    // Task 1.6. Distinguish four basis states
    // Input: two qubits (stored in an array) which are guaranteed to be
    //        in one of the four basis states (|00⟩, |01⟩, |10⟩ or |11⟩).
    // Output: 0 if qubits were in |00⟩ state,
    //         1 if they were in |01⟩ state,
    //         2 if they were in |10⟩ state,
    //         3 if they were in |11⟩ state.
    // In this task and the subsequent ones the order of qubit states
    // in task description matches the order of qubits in the array
    // (i.e., |10⟩ state corresponds to qs[0] in state |1⟩ and qs[1] in state |0⟩).
    // The state of the qubits at the end of the operation does not matter.
    operation BasisStateMeasurement (qs : Qubit[]) : Int {
        // ...
        mutable x = 0;
        if(M(qs[0]) == One){
            set x = x + 2;
        }
        if(M(qs[1]) == One){
            set x = x + 1;
        }
        return x;
    }
    
    
    // Task 1.7. Distinguish two basis states given by bit strings
    // Inputs:
    //      1) N qubits (stored in an array) which are guaranteed to be
    //         in one of the two basis states described by the given bit strings.
    //      2) two bit string represented as Bool[]s.
    // Output: 0 if qubits were in the basis state described by the first bit string,
    //         1 if they were in the basis state described by the second bit string.
    // Bit values false and true correspond to |0⟩ and |1⟩ states.
    // The state of the qubits at the end of the operation does not matter.
    // You are guaranteed that the both bit strings have the same length as the qubit array,
    // and that the bit strings will differ in at least one bit.
    // You can use exactly one measurement.
    // Example: for bit strings [false, true, false] and [false, false, true]
    //          return 0 corresponds to state |010⟩, and return 1 corresponds to state |001⟩.
    operation TwoBitstringsMeasurement (qs : Qubit[], bits1 : Bool[], bits2 : Bool[]) : Int {
        // ...
        let  n = Length(qs);
        for(i in 0..n-1){
            if(bits1[i] != bits2[i]){
                if((M(qs[i]) == One) == bits1[i]){
                    return 0;
                } else {
                    return 1;
                }
            }
        }
        return -1;
    }
    
    
    // Task 1.8. |0...0⟩ state or W state ?
    // Input: N qubits (stored in an array) which are guaranteed to be
    //        either in |0...0⟩ state
    //        or in W state (https://en.wikipedia.org/wiki/W_state).
    // Output: 0 if qubits were in |0...0⟩ state,
    //         1 if they were in W state.
    // The state of the qubits at the end of the operation does not matter.
    operation AllZerosOrWState (qs : Qubit[]) : Int {
        // ...
        let  n = Length(qs);
        for(i in 0..n-1){
            if(M(qs[i]) == One){ return 1;}
        }
        return 0;
    }
    
    
    // Task 1.9. GHZ state or W state ?
    // Input: N >= 2 qubits (stored in an array) which are guaranteed to be
    //        either in GHZ state (https://en.wikipedia.org/wiki/Greenberger%E2%80%93Horne%E2%80%93Zeilinger_state)
    //        or in W state (https://en.wikipedia.org/wiki/W_state).
    // Output: 0 if qubits were in GHZ state,
    //         1 if they were in W state.
    // The state of the qubits at the end of the operation does not matter.
    operation GHZOrWState (qs : Qubit[]) : Int {
        // ...
        let fs = M(qs[0]);
        let  n = Length(qs);
        for(i in 0..n-1){
            if(M(qs[i]) != fs){ return 1;}
        }
        return 0;
    }
    
    
    // Task 1.10. Distinguish four Bell states
    // Input: two qubits (stored in an array) which are guaranteed to be in one of the four Bell states:
    //         |Φ⁺⟩ = (|00⟩ + |11⟩) / sqrt(2)
    //         |Φ⁻⟩ = (|00⟩ - |11⟩) / sqrt(2)
    //         |Ψ⁺⟩ = (|01⟩ + |10⟩) / sqrt(2)
    //         |Ψ⁻⟩ = (|01⟩ - |10⟩) / sqrt(2)
    // Output: 0 if qubits were in |Φ⁺⟩ state,
    //         1 if they were in |Φ⁻⟩ state,
    //         2 if they were in |Ψ⁺⟩ state,
    //         3 if they were in |Ψ⁻⟩ state.
    // The state of the qubits at the end of the operation does not matter.
    operation BellPairTransform(qs: Qubit[]): Unit {
        body(...){
            H(qs[1]);
            CNOT(qs[1], qs[0]);
            // rev 0,1
            // 00 -> 00 + 10 -> 00 + 11 (0)
            // 10 -> 00 - 10 -> 00 - 11 (2)
            // 01 -> 01 + 11 -> 01 + 10 (1)
            // 11 -> 01 - 11 -> 01 - 10 (3)
            
        }
        adjoint invert;
    }
    operation BellState (qs : Qubit[]) : Int {
        // Hint: you need to use 2-qubit gates to solve this task
        
        Adjoint BellPairTransform(qs);
        return BasisStateMeasurement(qs);
    }
    
    
    // Task 1.11*. Distinguish four orthogonal 2-qubit states
    // Input: two qubits (stored in an array) which are guaranteed to be in one of the four orthogonal states:
    //         |S0⟩ = (|00⟩ + |01⟩ + |10⟩ + |11⟩) / 2
    //         |S1⟩ = (|00⟩ - |01⟩ + |10⟩ - |11⟩) / 2
    //         |S2⟩ = (|00⟩ + |01⟩ - |10⟩ - |11⟩) / 2
    //         |S3⟩ = (|00⟩ - |01⟩ - |10⟩ + |11⟩) / 2
    // Output: 0 if qubits were in |S0⟩ state,
    //         1 if they were in |S1⟩ state,
    //         2 if they were in |S2⟩ state,
    //         3 if they were in |S3⟩ state.
    // The state of the qubits at the end of the operation does not matter.
    operation Transform11(qs: Qubit[]): Unit {
        body(...){
            H(qs[1]);
            H(qs[0]);
            // rev 0,1
            // (0) 00 -> + + + +
            // (1) 01 -> + - + - 
            
        }
        adjoint invert;
    }
    operation TwoQubitState (qs : Qubit[]) : Int {
        Adjoint Transform11(qs);

        return BasisStateMeasurement(qs);
    }
    
    
    // Task 1.12**. Distinguish four orthogonal 2-qubit states, part two
    // Input: two qubits (stored in an array) which are guaranteed to be in one of the four orthogonal states:
    //         |S0⟩ = ( |00⟩ - |01⟩ - |10⟩ - |11⟩) / 2
    //         |S1⟩ = (-|00⟩ + |01⟩ - |10⟩ - |11⟩) / 2
    //         |S2⟩ = (-|00⟩ - |01⟩ + |10⟩ - |11⟩) / 2
    //         |S3⟩ = (-|00⟩ - |01⟩ - |10⟩ + |11⟩) / 2
    // Output: 0 if qubits were in |S0⟩ state,
    //         1 if they were in |S1⟩ state,
    //         2 if they were in |S2⟩ state,
    //         3 if they were in |S3⟩ state.
    // The state of the qubits at the end of the operation does not matter.
    operation Transform12(qs: Qubit[]): Unit {
        body(...){
            H(qs[0]);H(qs[1]);
            Z(qs[0]);
            X(qs[0]);
            Controlled Z([qs[0]],qs[1]);
            X(qs[0]);
        }
        adjoint invert;
    }
    operation TwoQubitStatePartTwo (qs : Qubit[]) : Int {
        // ...
        Adjoint Transform12(qs);
        mutable res = BasisStateMeasurement(qs);
        if(res == 2){ set res = 1;}
        elif(res == 1){ set res = 2;}
        return res;
    }
    
    
    //////////////////////////////////////////////////////////////////
    // Part II*. Discriminating Nonorthogonal States
    //////////////////////////////////////////////////////////////////
    
    // The solutions for tasks in this section are validated using the following method.
    // The solution is called on N input states, each of which is picked randomly,
    // with all possible input states equally likely to be generated.
    // The accuracy of state discrimination is estimated as an average of
    // discrimination correctness over all input states.
    
    // Task 2.1*. |0⟩ or |+⟩ ?
    //           (quantum hypothesis testing or state discrimination with minimum error)
    // Input: a qubit which is guaranteed to be in |0⟩ or |+⟩ state with equal probability.
    // Output: true if qubit was in |0⟩ state, or false if it was in |+⟩ state.
    // The state of the qubit at the end of the operation does not matter.
    // Note: in this task you have to get accuracy of at least 80%.
    operation IsQubitPlusOrZero (q : Qubit) : Bool {
        // ...
        Ry(PI()/4.0, q);
        let val = M(q);
        if(val == One) {return false;}
        else {return true;}
    }
    
    
    // Task 2.2**. |0⟩, |+⟩ or inconclusive?
    //             (unambiguous state discrimination)
    // Input: a qubit which is guaranteed to be in |0⟩ or |+⟩ state with equal probability.
    // Output: 0 if qubit was in |0⟩ state,
    //         1 if it was in |+⟩ state,
    //         -1 if you can't decide, i.e., an "inconclusive" result.
    // Your solution:
    //  - can never give 0 or 1 answer incorrectly (i.e., identify |0⟩ as 1 or |+⟩ as 0).
    //  - must give inconclusive (-1) answer at most 80% of the times.
    //  - must correctly identify |0⟩ state as 0 at least 10% of the times.
    //  - must correctly identify |+⟩ state as 1 at least 10% of the times.
    //
    // The state of the qubit at the end of the operation does not matter.
    // You are allowed to use ancilla qubit(s).
    operation IsQubitPlusZeroOrInconclusiveSimpleUSD (q : Qubit) : Int {
        let basis = RandomInt(2);
        if(basis == 1){ H(q);}
        let res = M(q);
        if(res == One && basis == 0){ return 1;}
        elif(res == One && basis == 1){return 0;}
        else {return -1;}
    }
    
}
