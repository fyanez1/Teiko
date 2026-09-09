-- Relational schema for immune cell-count data (SQLite).
--
-- Design summary (see README.md for the full rationale):
--   projects (1) --< subjects (1) --< samples (1) --< cell_counts >-- (1) cell_populations
--
-- * Metadata that is constant per subject (condition, age, sex, treatment,
--   response) lives once on `subjects`; per-collection facts (sample type,
--   time point) live on `samples`.
-- * Counts are stored in "long" form: one row per (sample, population).
--   Adding a new population is a data change, not a schema change.
-- * Views expose the analysis-ready shapes (per-sample relative frequencies
--   and a denormalised sample metadata join) so every consumer computes
--   percentages the same way.

PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS projects (
    project_id  TEXT PRIMARY KEY
);

CREATE TABLE IF NOT EXISTS subjects (
    subject_id  TEXT PRIMARY KEY,
    project_id  TEXT    NOT NULL REFERENCES projects(project_id),
    condition   TEXT    NOT NULL,                      -- melanoma | carcinoma | healthy | ...
    age         INTEGER CHECK (age IS NULL OR age >= 0),
    sex         TEXT    CHECK (sex IN ('M', 'F')),
    treatment   TEXT    NOT NULL,                      -- miraclib | phauximab | none | ...
    response    TEXT    CHECK (response IN ('yes', 'no'))  -- NULL when not applicable (e.g. healthy/untreated)
);

CREATE TABLE IF NOT EXISTS samples (
    sample_id                 TEXT PRIMARY KEY,
    subject_id                TEXT    NOT NULL REFERENCES subjects(subject_id),
    sample_type               TEXT    NOT NULL,        -- PBMC | WB | ...
    time_from_treatment_start INTEGER NOT NULL,        -- days
    UNIQUE (subject_id, sample_type, time_from_treatment_start)
);

CREATE TABLE IF NOT EXISTS cell_populations (
    population_id INTEGER PRIMARY KEY,
    name          TEXT    NOT NULL UNIQUE,             -- b_cell, cd8_t_cell, ...
    display_order INTEGER NOT NULL
);

CREATE TABLE IF NOT EXISTS cell_counts (
    sample_id     TEXT    NOT NULL REFERENCES samples(sample_id),
    population_id INTEGER NOT NULL REFERENCES cell_populations(population_id),
    count         INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample_id, population_id)
);

-- Indexes on the columns used for joins and filters.
CREATE INDEX IF NOT EXISTS idx_subjects_project   ON subjects(project_id);
CREATE INDEX IF NOT EXISTS idx_subjects_cohort    ON subjects(condition, treatment, response);
CREATE INDEX IF NOT EXISTS idx_samples_subject    ON samples(subject_id);
CREATE INDEX IF NOT EXISTS idx_samples_type_time  ON samples(sample_type, time_from_treatment_start);
CREATE INDEX IF NOT EXISTS idx_cell_counts_pop    ON cell_counts(population_id);

-- Total cells per sample (sum over all populations).
CREATE VIEW IF NOT EXISTS sample_totals AS
SELECT sample_id, SUM(count) AS total_count
FROM cell_counts
GROUP BY sample_id;

-- Part 2: relative frequency of every population in every sample.
CREATE VIEW IF NOT EXISTS sample_frequencies AS
SELECT
    c.sample_id                          AS sample,
    t.total_count                        AS total_count,
    p.name                               AS population,
    c.count                              AS count,
    100.0 * c.count / t.total_count      AS percentage,
    p.display_order                      AS display_order
FROM cell_counts AS c
JOIN cell_populations AS p ON p.population_id = c.population_id
JOIN sample_totals    AS t ON t.sample_id     = c.sample_id;

-- Denormalised sample metadata: one row per sample with its subject/project facts.
CREATE VIEW IF NOT EXISTS sample_metadata AS
SELECT
    s.sample_id,
    s.subject_id,
    sub.project_id,
    sub.condition,
    sub.age,
    sub.sex,
    sub.treatment,
    sub.response,
    s.sample_type,
    s.time_from_treatment_start
FROM samples  AS s
JOIN subjects AS sub ON sub.subject_id = s.subject_id;
