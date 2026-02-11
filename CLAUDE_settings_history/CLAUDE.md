# Default workflow (apply to every project)
- Do not implement code immediately.
- Always start with: 
	- (1) Planning phase (Codebase understanding if available, Architecture / design proposal), in exchange with the user (discussing design questions). 
	- (2) Implementation plan.
- All new findings and decisions are saved into markdown files.

## Phase 1 — Analysis (required if codebase existing)
- Produce an `ANALYSIS.md` (or a Markdown section in chat) that includes:
  - Problem restatement + constraints
  - Relevant files/modules to inspect
  - Risks/unknowns + assumptions
  - Proposed approach options (if applicable) with tradeoffs

## Phase 2 — Architecture (required)
- Propose the target architecture and interfaces before touching code:
  - Module boundaries
  - Data flow
  - Key types/contracts
  - Testing strategy
- Write an `ARCHITECTURE.md` with the proposed architecture.
- Discuss important design questions with the user and update the ARCHITECTURE.md.

## Phase 3 — Plan then implement
- Write a step-by-step implementation plan with checkpoints, store the implementation plan in `IMPLEMENTATION.md`.

## PHASE 4 - Implementation
- Implement the `IMPLEMENTATION.md` step by step. Between each step the user will review the code and trigger the next step.

