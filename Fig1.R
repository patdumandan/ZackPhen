# Load necessary libraries
library(ggplot2)
library(dplyr)

# Create data for the three scenarios
scenarios <- data.frame(
  Scenario = rep(c("Scenario 1: Matched timing", "Scenario 2: Prolonged activity duration, higher overlap",
                   "Scenario 3: Shortened activity duration: lower overlap"), each = 2),
  Population = rep(c("Consumer", "Resource"), 3),
  Start_Date = c(50, 50,   # Scenario 1: Both consumer and resource start at day 50
                 50, 40,   # Scenario 2: Resource starts earlier (day 40)
                 50, 60),  # Scenario 3: Resource starts later (day 60)
  End_Date = c(150, 150,   # Scenario 1: Both end at day 150
               150, 170,   # Scenario 2: Resource ends later (day 170)
               150, 130)   # Scenario 3: Resource ends earlier (day 130)
)

# Add a column for color-coding
scenarios$Color <- ifelse(scenarios$Population == "Consumer", "Consumer", "Resource")

# Create the ggplot
ggplot(scenarios, aes(x = Start_Date, xend = End_Date, y = Population,
                       color=Color)) +
  geom_segment(size = 4) +
  theme_minimal() +
  labs(title = "Phenology of Consumer and Resource Populations in Different Scenarios",
       x = "Day of Year",
       y = "Scenario",
       color = "Population") + facet_wrap(~Scenario, ncol=1, nrow=3)+
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_blank(),
        axis.title = element_text(size = 14),
        plot.title = element_text(hjust = 0.5, size = 16))
