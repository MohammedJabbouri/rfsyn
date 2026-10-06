#include "chains.h"
#include "chain.h"
#include <stdlib.h>
#include <string.h>

#define CHAINS_INITIAL_CAPACITY 8

typedef struct {
  char *name;
  transform_t *stages;
  size_t n_stages;
} chain_entry_t;

struct chains_registry {
  chain_entry_t *entries;
  size_t capacity;
  size_t count;
};

chains_registry_t *chains_registry_create(void) {
  
  chains_registry_t *reg = malloc(sizeof(chains_registry_t));
  if (!reg) return NULL;

  reg->entries = malloc(CHAINS_INITIAL_CAPACITY * sizeof(chain_entry_t));

  if (!reg->entries) {
    free(reg);
    return NULL;
  }

  reg->count = 0;
  reg->capacity = CHAINS_INITIAL_CAPACITY;
  return reg;
}

static void entry_release(chain_entry_t *e) {
  for (size_t i = 0; i < e->n_stages; i++) {
    if (e->stages[i].destroy) {
      e->stages[i].destroy(e->stages[i].ctx);
    }
  }
  free(e->stages);
  free(e->name);
}

static chain_entry_t *find_entry(const chains_registry_t *reg, const char *name) {
  for (size_t i = 0; i < reg->count; i++) {
      if (strcmp(reg->entries[i].name, name) == 0) {
        return &reg->entries[i];
      }
  }
  return NULL;
}

int chains_register(chains_registry_t *reg, const char *name,
const transform_t *stages, size_t n_stages) {
  if (!reg || !name || !stages || n_stages == 0) return -1;

// LIKELY TO FAIL! WILL FIX LATER
  
  size_t name_len = strlen(name) + 1;
  char *name_copy = malloc(name_len);
  if (!name_copy) return -1;
  memcpy(name_copy, name, name_len);

  transform_t *stages_copy = malloc(n_stages * sizeof(transform_t));

  if (!stages_copy) {
    free(name_copy);
    return -1;
  }
  
  memcpy(stages_copy, stages, n_stages * sizeof(transform_t));

  chain_entry_t *slot = find_entry(reg, name);
  
  if (slot) {
    entry_release(slot);
  } else {
    if (reg->count == reg->capacity) {
      size_t new_cap = reg->capacity * 2;
      chain_entry_t *grown = realloc(reg->entries, new_cap * sizeof(chain_entry_t));
      if (!grown) {
        free(name_copy);
        free(stages_copy);
        return -1;
      }
      reg->entries  = grown;
      reg->capacity = new_cap;
    }
    slot = &reg->entries[reg->count];
    reg->count++;
  }

  slot->name     = name_copy;
  slot->stages   = stages_copy;
  slot->n_stages = n_stages;

  return 0;
}

int chains_lookup(const chains_registry_t *reg, const char *name,
  const transform_t **stages_out, size_t *n_stages_out) {
    if (!reg || !name) return -1;
    chain_entry_t *e = find_entry(reg, name);
    if (!e) return -1;
    if (stages_out)   *stages_out   = e->stages;
    if (n_stages_out) *n_stages_out = e->n_stages;
    return 0;
}

int chains_run(const chains_registry_t *reg, const char *name, signal_t *sig) {
  const transform_t *stages;
  size_t n_stages;
  
  if (chains_lookup(reg, name, &stages, &n_stages) != 0) return -1;
  return chain_apply(stages, n_stages, sig);
}

void chains_registry_destroy(chains_registry_t *reg) {
  if (!reg) return;

  for (size_t i = 0; i < reg->count; i++) {
    entry_release(&reg->entries[i]);
  }
  free(reg->entries);
  free(reg);
}