from django.db import models

# Create your models here.
from django.contrib.auth.models import User
from django.utils.timezone import now

class InteractionResult(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)  # Links each result to a user
    input1 = models.TextField()
    input2 = models.TextField()
    result = models.CharField(max_length=20)  # Stores "Interacting" or "Non Interacting"
    timestamp = models.DateTimeField(auto_now_add=True)  # Auto add timestamp

    def __str__(self):
        return f"{self.user.username} - {self.result} ({self.timestamp})"


class InteractionResult_PP(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)  # Links each result to a user
    input1 = models.TextField()
    input2 = models.TextField()
    result = models.CharField(max_length=20)  # Stores "Interacting" or "Non Interacting"
    timestamp = models.DateTimeField(auto_now_add=True)  # Auto add timestamp

    def __str__(self):
        return f"{self.user.username} - {self.result} ({self.timestamp})"
    

class InteractionResult_PL(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)  # Links each result to a user
    input1 = models.TextField()
    input2 = models.TextField()
    result = models.CharField(max_length=20)  # Stores "Interacting" or "Non Interacting"
    timestamp = models.DateTimeField(auto_now_add=True)  # Auto add timestamp

    def __str__(self):
        return f"{self.user.username} - {self.result} ({self.timestamp})"