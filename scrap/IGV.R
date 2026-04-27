
dflig <- mutlig %>%
  dplyr::select(ref, ali, splcov, nreads, freq) %>%
  mutate(index=1:nrow(mutlig), .before=1)

ligl <- dflig %>%
  # pivot ref and ali
  pivot_longer(cols=c(ref, ali), names_to='type', values_to='seq') %>%
  # split sequences into single characters
  mutate(seqs=str_split(seq, '')) %>%
  unnest(seqs) %>%
  # position in sequence
  group_by(index, type) %>%
  mutate(pos=row_number()) %>%
  ungroup()

ligl <- ligl %>%
  filter(index%in%c(1, 2, 3))

ggseqs <- ggplot(ligl, aes(x=pos, y=type, label=seqs)) +
  geom_text(family="Courier", size = 4) +
  facet_wrap(~ index, ncol=1, scales="free_x") +
  theme_minimal(base_family = "Courier") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_blank(),
    axis.text.y = element_text(face = "bold"),
    panel.grid = element_blank()
  )
ggseqs
ggsave(here('250328_rab/plots/igvseqs.pdf'), ggseqs, width=500, height=500, units='mm')
