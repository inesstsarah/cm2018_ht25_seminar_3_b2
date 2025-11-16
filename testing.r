tumor_group$predicted_DV <- predict(m2_lme)

plot(tumor_group$DV, tumor_group$predicted_DV,
     xlab = "Actual DV", ylab = "Predicted DV",
     main = "Actual vs Predicted DV In Tumor Group",
     pch = 19, col = "darkgreen")
abline(0, 1, col = "red", lwd = 2)


tumor_group
# Extract random intercepts for ID
ranef_df <- ranef(m2_lme)$ID
ranef_df$ID <- rownames(ranef_df)
colnames(ranef_df) <- c("Intercept", "ID")
print(head(ranef_df))

# Extract random intercepts for ID
ranef_df <- ranef(m2_lme)$ID
ranef_df$ID <- rownames(ranef_df)
colnames(ranef_df) <- c("Intercept", "ID")


# Plot the random intercepts
ggplot(ranef_df, aes(x = reorder(ID, Intercept), y = Intercept)) +
  geom_point(color = "blue", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Random Intercepts by ID",
       x = "ID",
       y = "Intercept") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# Extract fixed effects
fixef_vals <- fixef(m2_lme)

# Extract random intercepts
ranef_vals <- ranef(m2_lme)$ID
ranef_df <- data.frame(ID = rownames(ranef_vals),
                       rand_intercept = ranef_vals[, "(Intercept)"])

# Ensure ID is character in both data frames
tumor_group$ID <- as.character(tumor_group$ID)
ranef_df$ID <- as.character(ranef_df$ID)

# Merge random effects into original data
tumor_group <- left_join(tumor_group, ranef_df, by = "ID")

# Compute group-specific intercepts and slopes
tumor_group$group_intercept <- fixef_vals[1] + tumor_group$rand_intercept
tumor_group$group_slope <- fixef_vals["Treatment"]

# Compute predicted values manually
tumor_group$predicted <- tumor_group$group_intercept + tumor_group$group_slope * tumor_group$Treatment

# Plot
ggplot(tumor_group, aes(x = Treatment, y = DV, color = ID)) +
  geom_point(alpha = 0.6) +
  geom_line(aes(y = predicted), size = 1) +
  geom_abline(aes(intercept = group_intercept, slope = group_slope, color = ID),
              data = distinct(tumor_group, ID, group_intercept, group_slope),
              linetype = "dashed", size = 0.8) +
  geom_abline(intercept = fixef_vals[1], slope = fixef_vals["Treatment"],
              color = "black", linetype = "solid", size = 1.2) +
  labs(title = "Group-Specific Lines with Fixed Effect Overlay",
       x = "Treatment", y = "DV") +
  theme_minimal()

# Compute predicted values manually
tumor_group$predicted <- tumor_group$group_intercept + tumor_group$group_slope * tumor_group$Treatment

# Plot: Time on x-axis
ggplot(tumor_group, aes(x = Time, y = DV, color = ID)) +
  geom_point(alpha = 0.6) +
  geom_line(aes(y = predicted), size = 1) +
  labs(title = "DV over Time by Group (with Predictions)",
       x = "Time (days)", y = "DV") +
  theme_minimal()
